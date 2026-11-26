"""Preflight input validation — run BEFORE the pipeline. Catches format problems
in seconds so a bad upload fails fast instead of 30 min into Foldseek. Portal calls
this on submit; CLI users can run it standalone:
    python workflow/scripts/validate_inputs.py config/config.yaml
Exit 0 = OK to run; exit 1 = blocking errors (printed). Warnings don't block.
"""
import os
import re
import sys

from id_utils import SAFE_ID_RE, extract_accession
from structure_io import list_structure_files, read_ca_metrics


def validate(cfg):
    errors, warns, info = [], [], {}
    strain = cfg.get("strain", {}).get("code", "")
    pdb_dir = cfg["input"]["pdb_dir"]
    custom_pattern = str(cfg.get("input", {}).get("accession_pattern", "") or "").strip() or None

    if custom_pattern:
        try:
            re.compile(custom_pattern)
        except re.error as exc:
            errors.append(f"input.accession_pattern is not a valid regex: {exc}")
            return errors, warns, info

    qc = cfg.get("qc", {})
    try:
        min_length = int(qc.get("min_length", 50))
        max_length = int(qc.get("max_length", 1000))
        min_plddt = float(qc.get("min_plddt", 50))
    except (TypeError, ValueError):
        errors.append("qc.min_length, qc.max_length, and qc.min_plddt must be numeric")
        return errors, warns, info
    if min_length < 1 or max_length < 1:
        errors.append("QC protein-length limits must be at least 1 aa")
    if min_length > max_length:
        errors.append("qc.min_length cannot exceed qc.max_length")
    if not 0 <= min_plddt <= 100:
        errors.append("qc.min_plddt must be between 0 and 100")
    if min_length < 40:
        warns.append(
            "qc.min_length is below 40 aa; very short chains can have less "
            "discriminating TM-scores"
        )

    # ---- structures (required: .pdb, .cif, .mmcif, optional .gz) ----
    structures = list_structure_files(pdb_dir)
    info["n_structures"] = len(structures)
    if not structures:
        errors.append(f"no .pdb/.cif/.mmcif files in {pdb_dir}")
    bad_name, no_ca, lens, dup_accs = [], [], [], []
    seen_acc = {}
    accs = set()
    mmcif_count = 0
    for p in structures[:100000]:
        b = p.name
        if b.lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
            mmcif_count += 1
        acc = extract_accession(
            b,
            strain=strain,
            custom_pattern=custom_pattern,
            strip_strain_prefix=True,
        )
        if not acc or not SAFE_ID_RE.match(acc):
            bad_name.append(b)
        else:
            if acc in seen_acc:
                dup_accs.append(f"{acc} ({seen_acc[acc]}, {b})")
            else:
                seen_acc[acc] = b
                accs.add(acc)
        try:
            ca, _ = read_ca_metrics(p)
        except Exception:
            ca = 0
        if ca == 0:
            no_ca.append(b)
        else:
            lens.append(ca)
    if mmcif_count:
        info["n_mmcif_structures"] = mmcif_count
    if bad_name:
        errors.append(
            f"{len(bad_name)} structure filenames could not yield a valid protein ID: "
            f"e.g. {bad_name[:3]}"
        )
    if dup_accs:
        errors.append(
            f"{len(dup_accs)} duplicate protein accession(s) detected across structures: "
            f"e.g. {dup_accs[:3]}"
        )
    if no_ca:
        errors.append(f"{len(no_ca)} structures have no CA atoms (empty/corrupt): e.g. {no_ca[:3]}")
    if lens:
        under = [l for l in lens if l < min_length]
        over = [l for l in lens if l > max_length]
        if under:
            warns.append(
                f"{len(under)} structures are shorter than min_length "
                f"{min_length} aa (will fail QC)"
            )
        if over:
            warns.append(
                f"{len(over)} structures exceed max_length "
                f"{max_length} aa (will fail QC)"
            )
        info["length_range"] = [min(lens), max(lens)]
    if strain and structures and not any(p.name.startswith(strain + "_") for p in structures[:50]):
        warns.append(
            f"no structure filename starts with strain prefix '{strain}_' "
            "— filenames will be resolved by extracted protein accession"
        )

    # ---- sequences (optional) ----
    seqs_fa = cfg["input"].get("seqs_fasta", "")
    if seqs_fa and os.path.exists(seqs_fa):
        seq_accs = set()
        n_seq = 0
        with open(seqs_fa, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    n_seq += 1
                    acc = extract_accession(
                        line[1:],
                        strain=strain,
                        custom_pattern=custom_pattern,
                        strip_strain_prefix=False,
                    )
                    if acc:
                        seq_accs.add(acc)
        info["n_sequences"] = n_seq
        missing = accs - seq_accs
        if missing:
            warns.append(
                f"{len(missing)} structures have no matching sequence "
                f"(ESM/annotation will derive from structure): e.g. {sorted(missing)[:3]}"
            )
    elif seqs_fa:
        warns.append(f"seqs_fasta '{seqs_fa}' not found — sequences derived from structures")

    # ---- RNAseq (optional, 2-sheet standard format) ----
    rna = cfg["input"].get("rnaseq_xlsx", "")
    if rna and os.path.exists(rna):
        try:
            import openpyxl

            wb = openpyxl.load_workbook(rna, read_only=True)
            sheets = wb.sheetnames
            info["rnaseq_sheets"] = sheets
            need = {"id_mapping", "expression"}
            miss = need - set(sheets)
            if miss:
                errors.append(
                    f"RNAseq xlsx missing required sheet(s): {sorted(miss)} "
                    "(need id_mapping + expression) — see RNAseq_template"
                )
            else:
                idm = wb["id_mapping"]
                hdr = [c.value for c in next(idm.iter_rows(max_row=1))]
                if not (hdr[:2] == ["protein_accession", "gene_id"]):
                    errors.append(
                        "id_mapping sheet headers must be "
                        f"[protein_accession, gene_id], got {hdr[:2]}"
                    )
                else:
                    mapped = sum(
                        1 for r in idm.iter_rows(min_row=2, values_only=True) if r and r[0]
                    )
                    info["rnaseq_mapped_accessions"] = mapped
                    exp = wb["expression"]
                    ehdr = [c.value for c in next(exp.iter_rows(max_row=1))]
                    if not ehdr or ehdr[0] != "gene_id":
                        errors.append("expression sheet first column must be 'gene_id'")
                    ncond = len([h for h in ehdr[1:] if h])
                    if ncond < 2:
                        warns.append(
                            f"expression sheet has only {ncond} sample column(s); "
                            "need >=2 for condition grouping"
                        )
        except Exception as e:
            errors.append(f"cannot read RNAseq xlsx ({type(e).__name__}: {e})")
    elif rna:
        errors.append(f"rnaseq_xlsx '{rna}' set in config but file not found")
    else:
        info["rnaseq"] = "none — RNAseq step will be skipped"

    return errors, warns, info


if __name__ == "__main__":
    import yaml

    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
    cfg = yaml.safe_load(open(cfg_path))
    errors, warns, info = validate(cfg)
    print("=== SUSS input preflight ===")
    for k, v in info.items():
        print(f"  {k}: {v}")
    for w in warns:
        print(f"  ⚠ WARN: {w}")
    for e in errors:
        print(f"  ✗ ERROR: {e}")
    if errors:
        print(f"\nBLOCKED — {len(errors)} error(s). Fix before running.")
        sys.exit(1)
    print(f"\nOK to run{' (' + str(len(warns)) + ' warning(s))' if warns else ''}.")
    sys.exit(0)

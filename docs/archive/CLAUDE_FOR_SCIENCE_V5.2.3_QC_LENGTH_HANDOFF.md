# SUSS Protein Atlas v5.2.3 QC-length handoff

## Why this release exists

The Botrytis cinerea run `20260811-014819-bc-5fa1ba57` contained 17 structures
shorter than the former fixed 50-aa portal default. Users could not change the
QC range, and domain Foldseek searched all uploaded structures rather than the
QC-passing subset. A failed-QC protein could therefore enter a D family but be
absent from the canonical sequence FASTA, terminating `domain_sequences`.

## Changes

- The portal exposes minimum protein length, maximum protein length, and minimum
  mean pLDDT. Defaults remain 50 aa, 1000 aa, and 50.
- Values are validated server-side, displayed on the result page, stored in the
  job manifest, and written to `config.yaml` under `qc`.
- A new `qc_structures` rule creates `results/qc_pdb/` from `qc.csv`. Both
  full-length and domain Foldseek now search exactly this directory.
- Domain sequence extraction warns and skips a missing parent defensively. This
  is a fallback; the shared QC structure directory is the primary fix.
- Preflight reports structures below and above the selected range and warns when
  `min_length < 40`, because TM-score becomes less discriminating for very short
  chains.

## Scientific interpretation

The length setting does not trim FASTA sequences. It decides which mature
structures enter every downstream analysis and is measured from CA residues.
Short proteins can be biologically meaningful and are no longer blocked by an
uneditable deployment default. Researchers should still report the cutoff and
inspect very short-chain structural links carefully.

## Claude acceptance checks

1. Set minimum length to 20 aa in the portal and confirm the generated config
   contains `qc.min_length: 20`.
2. Submit a mixture containing one QC-pass and one QC-fail structure. Confirm
   `results/qc_pdb_manifest.csv` contains only the passing filename.
3. Confirm both Foldseek commands use `results/qc_pdb`.
4. Confirm an inverted range, such as minimum 100 and maximum 20, is rejected
   before a job is queued.
5. Review families containing chains below 40 aa with the short-chain warning in
   mind; this release changes inclusion policy, not the TM-score definition.

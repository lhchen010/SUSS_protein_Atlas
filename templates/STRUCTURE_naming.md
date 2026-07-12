# Structure File Naming Rules

Uploaded AF2 structure PDB filenames must be:

    <organism_code>_<accession>.pdb

Example: `cor_TDZ24916.1.pdb` (cor = C. orbiculare, TDZ24916.1 = GenBank accession)

## Rules
- **organism_code** is entered once at the upload interface, or written directly into the filename prefix. When merging multiple strains, avoid accession number collisions.
- **accession** must be recognized by `[A-Z]{2,3}\d{4,}\.\d+` (GenBank protein accession format).
- Structure must be **mature** (signal peptide removed, pre-processed before AF2 submission), pLDDT>50, length 50–1000 aa.
- One protein per .pdb file; may be packaged as zip for upload.

## Preflight checks run automatically before upload
- Filename contains valid accession
- PDB has CA atoms (empty/corrupted files rejected)
- Length does not exceed limit
- (if sequence provided) sequence accession matches structure
- (if RNAseq provided) xlsx is standard two-sheet format, gene_id matches

Errors block execution; warnings proceed anyway.

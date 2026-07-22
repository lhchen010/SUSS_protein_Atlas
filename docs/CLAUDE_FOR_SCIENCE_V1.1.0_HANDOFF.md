# SUSS Protein Atlas v1.1.0: Claude for Science handoff

## Scope

This release adds client-side family-network search and a mobile network layout. It does not change
clustering, family membership, annotation calls, TM-score calculations, SUSS classification, or
pipeline outputs.

## Search behaviour

The network toolbar now performs case-insensitive, multi-term search across:

- family ID and member accession/gene name;
- family consensus annotation and every member's Pfam/InterPro annotation;
- PDB and AFDB-SwissProt Foldseek hits and protein names;
- EffectorP result;
- DeepTMHMM transmembrane-region count;
- novelty status, structural TM, SUSS, identity, pLDDT, and length metrics.

Plain terms search every field. Precise filters accept `gene:`, `acc:`, `annotation:`, `effectorp:`,
`tmr:`, `deeptmhmm:`, `structtm:`, `family:`, `novel:`, and `suss:`. Multiple terms use AND
semantics.

Matching nodes keep their original scientific fill and receive an orange outline and shadow.
Non-matching nodes and edges are muted but remain visible for network context. Enter opens a unique
match or fits multiple matches; Escape and the clear button restore the full network.

## Payload change

Per-member network annotation objects now additionally retain `gene`, complete InterPro entries,
and the AFDB accession. Existing payload fields and workbook exports are unchanged.

## Responsive layout

At widths below 700 px, the network and family panel are stacked vertically. This removes the old
520 px horizontal overflow while preserving the network as the first interactive view.

## Staging validation

The builder regenerated the six-family, 100-protein production dataset from job
`20260714-100412-cor-436dde93`. Browser acceptance results:

| Query | Expected result | Observed result |
|---|---|---|
| `gene:TDZ13877.1` | F1 | 1 cluster; Enter opened F1 |
| `annotation:Peroxidase` | F5 | 1 cluster |
| `effectorp:non-effector` | F0, F1, F2, F5 | 4 clusters |
| `tmr:1` | F1, F2 | 2 clusters |

Desktop highlighting was visually inspected. A 390 x 844 px browser check reported a 390 px body,
390 px network, no horizontal overflow, and the correct two-cluster TMR result.

## Automated validation

- Six runtime/portal unit tests passed.
- Seven atlas tests passed, including the new search UI and field contract.
- Python compile and whitespace checks passed.
- GitHub Actions must run the complete 13-test suite before merge.

## Claude acceptance checks

1. Search for a member accession and confirm only its family remains emphasized.
2. Search `annotation:Peroxidase` and confirm F5 is the only match.
3. Search `effectorp:non-effector` and confirm four clusters match.
4. Search `tmr:1` and confirm F1 and F2 match.
5. Press Enter on the accession query and confirm F1 opens in the side panel.
6. Clear the query and confirm all six nodes return to their original colors.
7. Repeat an annotation search on a narrow browser window and confirm the network remains visible.

Release version: `1.1.0`.

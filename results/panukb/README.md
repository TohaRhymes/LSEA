# Pan-UKB Aggregated Results

Lightweight aggregated results from LSEA analysis on Pan-UK Biobank data.

## Files

| File | Gene Set Category | Description |
|------|------------------|-------------|
| `C2__associated_phenos.tsv` | MSigDB C2 (KEGG) | Enriched phenotype-pathway associations |
| `GTE__associated_phenos.tsv` | GTEx tissues | Enriched phenotype-tissue associations |
| `BCM__associated_phenos.tsv` | Blood cell markers | Enriched phenotype-cell type associations |
| `GO:BP__associated_phenos.tsv` | GO Biological Process | Enriched phenotype-GO BP associations |
| `GO:CC__associated_phenos.tsv` | GO Cellular Component | Enriched phenotype-GO CC associations |
| `GO:MF__associated_phenos.tsv` | GO Molecular Function | Enriched phenotype-GO MF associations |

## How to obtain

Copy from the server path `${PANUKB_LSEA_DIR}/aggregated_checks/` (see `.env`).

## Context

These TSV files summarize which phenotypes show significant enrichment
for each gene set category (q-value < 0.05, >= 3 overlapping loci).

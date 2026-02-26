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

## Format

Columns: `gene_set`, `associated_phenotypes`, `phenos`

- `gene_set`: gene set / pathway name
- `associated_phenotypes`: number of phenotypes with significant enrichment
- `phenos`: comma-separated list of phenotype names

## How to obtain

Data files are not tracked in git. Copy from the server:
```bash
scp ${PANUKB_LSEA_DIR}/aggregated_checks/*__associated_phenos.tsv results/panukb/
```

## Context

These TSV files summarize which phenotypes show significant enrichment
for each gene set category (q-value < 0.05, >= 3 overlapping loci).

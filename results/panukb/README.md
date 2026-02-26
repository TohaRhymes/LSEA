# Pan-UKB Aggregated Results

Aggregated results from LSEA enrichment analysis (Experiment III) on 150 Pan-UK Biobank
phenotypes (139 unique phenotype codes) from the maximal independent set (MIS).

## Files

| File | Gene Set Category | Sets | Description |
|------|------------------|------|-------------|
| `C2__associated_phenos.tsv` | MSigDB C2 (all curated) | 7,233 | Enriched phenotype-pathway associations |
| `GTE__associated_phenos.tsv` | GTEx v8 tissues | 45 | Enriched phenotype-tissue associations |
| `BCM__associated_phenos.tsv` | Blood cell markers | 12 | Enriched phenotype-cell type associations |
| `GO:BP__associated_phenos.tsv` | GO Biological Process | 7,608 | Enriched phenotype-GO BP associations |
| `GO:CC__associated_phenos.tsv` | GO Cellular Component | 1,026 | Enriched phenotype-GO CC associations |
| `GO:MF__associated_phenos.tsv` | GO Molecular Function | 1,820 | Enriched phenotype-GO MF associations |

## Format

Columns: `gene_set`, `associated_phenotypes`, `phenos`

- `gene_set`: gene set / pathway name
- `associated_phenotypes`: number of phenotypes with significant enrichment
- `phenos`: comma-separated list of phenotype names

Significance criteria: q-value < 0.05, at least 3 overlapping loci.

## How to obtain

Data files are not tracked in git. Copy from the server:
```bash
scp ${PANUKB_LSEA_DIR}/aggregated_checks/*__associated_phenos.tsv results/panukb/
```

See `experiments/README.md` for full experiment details and `experiments/DATA_PREPARATION.md`
for data provenance.

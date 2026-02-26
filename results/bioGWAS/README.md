# bioGWAS Aggregated Results

Aggregated results from LSEA analysis on bioGWAS simulated data (Experiments I-II).
Includes comparison with MAGMA (linreg, mean, top) and PASCAL.

## Files

| File | Description |
|------|-------------|
| `TPR_to_draw_LSEA.csv` | True Positive Rate — continuous traits (50 iter x 4 pathway sizes) |
| `FPR_to_draw_LSEA.csv` | False Positive Rate — continuous traits |
| `binTPR_to_draw_LSEA.csv` | True Positive Rate — binary traits (50 iter x 4 pathway sizes) |
| `binFPR_to_draw_LSEA.csv` | False Positive Rate — binary traits |

## Format

Columns: `model`, `path`, `score`, `min`, `max`

- `model`: enrichment method — `LSEA`, `linreg` (MAGMA), `mean` (MAGMA), `top` (MAGMA), `PASCAL`
- `path`: pathway size — `path_small` (17 genes), `path_medium` (69), `path_big` (199), `path_random` (control)
- `score`: rate (TPR or FPR)
- `min`, `max`: 95% confidence interval bounds

## How to obtain

Data files are not tracked in git. Copy from the server:
```bash
scp ${LSEA_TEST_DIR}/aggregated_data/*_to_draw_LSEA.csv results/bioGWAS/
```

## Context

These CSV files are produced by post-processing notebooks that aggregate
`annotation_stats_*.tsv` files from individual LSEA runs into TPR/FPR
metrics. Simulations use N=10,000 individuals, K=30 causal SNPs (k=15 from
target KEGG pathway + 15 random). Gene sets: MSigDB C2 KEGG (186 pathways).

See `experiments/README.md` for full experiment details and `experiments/DATA_PREPARATION.md`
for data provenance.

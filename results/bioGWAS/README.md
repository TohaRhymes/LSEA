# bioGWAS Aggregated Results

Lightweight aggregated results from LSEA analysis on bioGWAS simulated data.

## Files

| File | Description |
|------|-------------|
| `TPR_to_draw_LSEA.csv` | True Positive Rate — continuous traits (50 iter x 4 pathway sizes) |
| `FPR_to_draw_LSEA.csv` | False Positive Rate — continuous traits |
| `binTPR_to_draw_LSEA.csv` | True Positive Rate — binary traits |
| `binFPR_to_draw_LSEA.csv` | False Positive Rate — binary traits |

## Format

Columns: `model`, `path`, `score`, `min`, `max`

- `model`: LSEA, linreg, mean, top, PASCAL
- `path`: pathway size (path_small, path_medium, path_big, path_random)
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
metrics for comparison with MAGMA and PASCAL.

# bioGWAS Aggregated Results

Lightweight aggregated results from LSEA analysis on bioGWAS simulated data.

## Files

| File | Description |
|------|-------------|
| `TPR_to_draw_LSEA.csv` | True Positive Rate — continuous traits (50 iter x 4 pathway sizes) |
| `FPR_to_draw_LSEA.csv` | False Positive Rate — continuous traits |
| `binTPR_to_draw_LSEA.csv` | True Positive Rate — binary traits |
| `binFPR_to_draw_LSEA.csv` | False Positive Rate — binary traits |

## How to obtain

Copy from the server path defined as `LSEA_TEST_DIR` in your `.env` file.

## Context

These CSV files are produced by post-processing notebooks that aggregate
`annotation_stats_*.tsv` files from individual LSEA runs into TPR/FPR
metrics for comparison with MAGMA and PASCAL.

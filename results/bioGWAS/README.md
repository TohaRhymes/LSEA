# bioGWAS Aggregated Results

Lightweight aggregated results from LSEA analysis on bioGWAS simulated data.

## Files

Copy these from the server (`/media/DATA/gwasim/round2/lsea_test/`):

| File | Description |
|------|-------------|
| `TPR_to_draw_LSEA.csv` | True Positive Rate — continuous traits (50 iter x 4 pathway sizes) |
| `FPR_to_draw_LSEA.csv` | False Positive Rate — continuous traits |
| `binTPR_to_draw_LSEA.csv` | True Positive Rate — binary traits |
| `binFPR_to_draw_LSEA.csv` | False Positive Rate — binary traits |

## How to obtain

```bash
# From the server:
SRC=/media/DATA/gwasim/round2/lsea_test
scp server:${SRC}/TPR_to_draw_LSEA.csv .
scp server:${SRC}/FPR_to_draw_LSEA.csv .
scp server:${SRC}/binTPR_to_draw_LSEA.csv .
scp server:${SRC}/binFPR_to_draw_LSEA.csv .
```

## Context

These CSV files are produced by post-processing notebooks that aggregate
`annotation_stats_*.tsv` files from individual LSEA runs into TPR/FPR
metrics for comparison with MAGMA and PASCAL.

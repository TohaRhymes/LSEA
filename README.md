# LSEA v0.2.4

*This is a new version of the tool that is now being extensively benchmarked. Please note that this is a development repository. Please use the contents of this repository with caution.*

---

## Prerequisites

- **PLINK** v1.9+
- **BEDTools** v2.30+
- **Python** (tested on 3.11): see `requirements.txt`

Install Python dependencies:
```bash
pip install -r requirements.txt
```

---

## Overview

LSEA (Locus Set Enrichment Analysis) is a tool for gene set enrichment analysis using GWAS summary statistics, accounting for linkage disequilibrium (LD) and locus structure.  
It allows you to:
- **Generate a universe** of genomic intervals and feature mappings (once per reference/annotation).
- **Run enrichment analysis** for new GWAS datasets using the precomputed universe.

---

## Usage

### 1. Universe Generation

You need to generate a "universe" file (JSON) describing all possible intervals and their mapping to features (genes, peaks, etc).  
There are two main modes:

#### **A. Using a single BED + GMT file**
```bash
python3 universe_generator.py \
    -v ./data/positions.tsv \
    -ft ./data/gencode_formatted.bed ./data/c2.all.v7.0.symbols.gmt \
    -i 100000 \
    -o ./data/universe.json
```
- `-v`/`--variants`: TSV file with at least three columns: chromosome, position, and variant ID (must be unique).
- `-ft`/`--features`: BED file with all features and GMT file with gene set definitions.
- `-i`/`--interval`: Window size (bp) around each variant (default: 500000).
- `-o`/`--out_json`: Output JSON file.

#### **B. Using a directory of per-set BED files**
```bash
python3 universe_generator.py \
    -v ./data/positions.tsv \
    -ffdir ./data/bed_sets/ \
    -i 100000 \
    -o ./data/universe_dir.json
```
- `-ffdir`/`--feature_files_dir`: Directory with one BED file per feature set.

**Additional options:**
- `-vc`/`--variants_colnames`: Specify column names for the variants file (default: chr pos id).
- `-tmp`/`--keep_temp`: Keep temporary files for debugging.

**Full documentation and all options are available here:**
```
python3 universe_generator.py --help
```

---

### 2. Running LSEA Enrichment

Once you have a universe file, you can run LSEA on your GWAS summary statistics:

```bash
python3 LSEA_2.4.py \
    --input ./data/in.tsv \
    --universe ./data/universe.json \
    --plink_dir <plink_folder_path> \
    --bfile <plink_prefix> \
    --out ./results_lsea
```

**Key options:**
- `--input`/`-i`: Input GWAS summary statistics (TSV, must have header with chromosome, position, ID, p-value).
- `--universe`/`-u`: One or more universe JSON files (space-separated).
- `--plink_dir`/`-pd`: Path to PLINK executable.
- `--bfile`/`-bf`: Prefix for PLINK binary genotype files (bed/bim/fam).
- `--out`/`-o`: Output directory (default: lsea_result).
- `--column_names`/`-cols`: Specify column names for input TSV (default: chr pos id p).
- `--clump_p1`/`-c_p`: List of p-value cutoffs for clumping (default: 1e-5 5e-8).
- `--clump_p2`, `--clump_r2`, `--clump_kb`: Advanced PLINK clumping parameters.
- `--qval_threshold`/`-qt`: Q-value threshold for reporting (default: 0.05; set to 1.0 combined with `--print_all` to see everything).
- `--interval_count_threshold`/`-ict`: Minimum number of intervals (loci) a gene set must overlap to be reported (default: 3).
- `--print_all`/`-a`: Output all results, not just significant ones.
- `--keep_temp`/`-tmp`: If set, keep intermediate files (mostly for debugging). By default, intermediate files are deleted.

**Full documentation and all options are available here:**
```
python3 LSEA_2.4.py --help
```

---

## Output

- All results are written to the output directory (default: `lsea_result`).
- For each universe and p-value cutoff (by default, LSEA uses two p-value cutoffs (`1e-05, 5e-08`); if you want to test more cutoffs, specify them using the `--clump_p1` option separated by space), you get:
  - `*_result_<cutoff>.tsv`: Enrichment results for each gene set.
  - `annotation_stats_<universe>.tsv`: Summary statistics for each run.
- Log messages and warnings are printed to the console for transparency and debugging.

---

## Error Handling & Logging

- The tool logs warnings for any malformed or missing lines/files during parsing.
- All critical errors (missing files, wrong formats) are logged and cause the program to exit with an error message.
- Temporary files are cleaned up unless `--keep_temp` is specified (for LSEA_2.4.py and universe_generator.py).

---

## Example Workflow

1. **Generate universe:**
   ```bash
   python3 universe_generator.py -v ./data/positions.tsv -ft ./data/gencode_formatted.bed ./data/c2.all.v7.0.symbols.gmt -i 100000 -o ./data/universe.json
   ```

2. **Run LSEA:**
   ```bash
   python3 LSEA_2.4.py --input ./data/in.tsv --universe ./data/universe.json --plink_dir /path/to/plink --bfile /path/to/genotypes --out ./results_lsea
   ```

---

## Troubleshooting

- If you see warnings about skipped lines, check your input files for formatting issues.
- If you get errors about missing columns or files, verify your file paths and headers.
- For PLINK or BEDTools errors, ensure they are installed and available in your PATH.

---

## Contact

For questions, bug reports, or contributions, please open an issue or contact the maintainers.

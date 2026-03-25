# Data Preparation for LSEA Experiments

This document describes how the input data for LSEA experiments was prepared.

## 1. Common Reference Data

### 1.1 Genotype Data (LD Reference Panel)

1000 Genomes Project, EUR population, merged into PLINK format:

```
merged_1000genomes_eur.{bed,bim,fam}
```

Used for:
- bioGWAS simulations (genotype source)
- PLINK LD-based clumping in both bioGWAS and Pan-UKB experiments

### 1.2 Gene Annotations

GENCODE v37 (`gencode.v37.annotation.gtf`) — used to map genes to genomic coordinates.

Converted to BED format (gene-level) with:
```bash
awk -F'\t' 'NR>=6 && $3=="gene"' gencode.v37.annotation.gtf | \
    awk '{sub(/^chr/, "", $1); match($0, /gene_name "([^"]+)"/, arr);
          geneName=arr[1]; print $1"\t"$4"\t"$5"\t"geneName}' > anno.bed
```

### 1.3 Gene Set Definitions (GMT Format)

| Short name | Source | File | Used in |
|-----------|--------|------|---------|
| `c2` (bioGWAS) | MSigDB C2 KEGG | `c2.cp.kegg.v2023.1.Hs.symbols.gmt` | bioGWAS |
| `c2` (Pan-UKB) | MSigDB C2 (all) | `c2.all.v2023.2.Hs.symbols.gmt` | Pan-UKB |
| `gte` | GTEx v8 | `GTEx8_formatted.gmt` | Pan-UKB |
| `bcm` | Literature | `blood_cell_markers.gmt` | Pan-UKB |
| `go_bp` | MSigDB C5 | `c5.go.bp.v2024.1.Hs.symbols.gmt` | Pan-UKB |
| `go_cc` | MSigDB C5 | `c5.go.cc.v2024.1.Hs.symbols.gmt` | Pan-UKB |
| `go_mf` | MSigDB C5 | `c5.go.mf.v2024.1.Hs.symbols.gmt` | Pan-UKB |

---

## 2. bioGWAS Simulated Data

Simulations were generated using [bioGWAS](https://github.com/TohaRhymes/GWAS_simulator)
(run via Docker).

### 2.1 Pathway Selection

Four pathway sizes were chosen from MSigDB C2 KEGG:

| Size | Pathway | Genes |
|------|---------|-------|
| `small` | KEGG_STEROID_BIOSYNTHESIS | 17 |
| `medium` | KEGG_PPAR_SIGNALING_PATHWAY | 69 |
| `big` | KEGG_FOCAL_ADHESION | 199 |
| `random` | ALL_GENES (control) | all |

Pathways were selected in `0_path_pick.ipynb`.

`ALL_GENES` was created with:
```bash
echo -e "ALL_GENES\thttps://github.com/TohaRhymes/GWAS_simulator\t$(awk '{print $1}' \
    gencode.v37.annotation.gtf.loc | sort | uniq | paste -s -d '\t')" > all.gmt
```

### 2.2 Simulation Parameters

All simulations use N=10,000 individuals, K=30 total causal SNPs.

| Experiment | Script | Iterations | k | Output dir | Pattern |
|-----------|--------|-----------|---|-----------|---------|
| Continuous (0-29) | `1_iterate_pathways.py` | 30 | 15 (30 for random) | `in_data/` | `test10000_path_{size}_{i}` |
| Continuous (30-49) | `extra_1_iterate_pathways.py` | 20 | 15 (30 for random) | `extra_in_data/` | `test10000_path_{size}_{i}` |
| Binary (0-49) | `binary_1_iterate_pathways.py` | 50 | 15 (30 for random) | `binary_in_data/` | `bin10000_path_{size}_{i}` |
| Sensitivity | `k_1_iterate_pathways.py` | 3 per k | 1..15 | `extra_in_data/` | `k10000_path_{size}_k{k}_{i}` |

**Simulation model parameters:**

| Parameter | Iter 0-29 | Iter 30-49 / Binary / Sensitivity |
|-----------|----------|-----------------------------------|
| `m_beta` | 0.05 | 0.5 |
| `sd_beta` | 0.001 | 0.01 |
| `gen_var` | 0.1 | 0.5 |
| `alpha` | (default) | 0.5 |
| `theta` | (default) | 1 |
| Genotypes | Simulated fresh | Pre-existing (`--skip_simulation --bfile_in_flag`) |

### 2.3 Output Format

Each simulation produces a GWAS summary statistics file:
```
{pattern}_{causal_id}_gwas.tsv
```
Columns include: `chr`, `pos`, `rsid`, `pval`, plus effect sizes and other statistics.

### 2.4 Universe Preparation

Variant positions are extracted from any GWAS file (they all share the same genotype panel):

```bash
awk '{print $1"\t"$3"\t"$2}' some_gwas.tsv > variants.tsv
```

Then the universe is built once using `universe_generator.py` (see `scripts/1_1_create_universe.sh`).

---

## 3. Pan-UK Biobank Data

### 3.1 Phenotype Selection

1. Download phenotype manifests from [Pan-UKB](https://pan.ukbb.broadinstitute.org/):
   - `phenotype_manifest.tsv` (phenotype metadata)
   - `h2_manifest.tsv` (heritability estimates)

2. Filter phenotypes (script: `_select_panukb_h2.py`):
   - Population: EUR
   - Sex: both_sexes
   - Condition: `in_max_independent_set == True`
   - Select phenotype with highest `h2_z` per unique `phenocode`

3. Download summary statistics:
   ```bash
   python3 download_panukb_summstats.py phenotype_manifest.tsv
   ```
   Downloads `.bgz` files from AWS links.

### 3.2 Decompress and Extract Columns

Script: `mv_to_dir_and_decompress_bgz.sh`

```bash
bgzip -cd file.bgz | awk 'BEGIN{FS=OFS="\t"}
    NR==1 { for (i=1; i<=NF; i++) f[$i]=i }
    { print $f["chr"], $f["pos"], $f["ref"], $f["alt"],
            $f["beta_EUR"], $f["se_EUR"], $f["neglog10_pval_EUR"] }
' > ukb_summstats/file.tsv.tsv
```

Result: 150 files in `ukb_summstats/*.tsv.tsv` with columns:
`chr`, `pos`, `ref`, `alt`, `beta_EUR`, `se_EUR`, `neglog10_pval_EUR`

### 3.3 Normalize for LSEA

Script: `preproc_tsv.py` (also available in `tests/preprocessing/preproc_tsv.py`)

Transformations:
1. Convert `neglog10_pval_EUR` to `pval` (`= 10^(-neglog10_pval_EUR)`)
2. Create `rsid` column as `chr:pos:ref:alt`
3. Filter to variants present in the 1000G EUR BIM file

```bash
# Single file:
python3 preproc_tsv.py INPUT.tsv.tsv OUTPUT.norm.tsv /path/to/merged_1000genomes_eur.bim

# Batch processing:
bash run_preproc.sh
```

Result: 150 files in `ukb_summstats/*.norm.tsv` with columns:
`chr`, `pos`, `ref`, `alt`, `beta_EUR`, `se_EUR`, `pval`, `rsid`

### 3.4 Universe Preparation

Variant positions are extracted from any Pan-UKB GWAS file:
```bash
awk 'BEGIN {OFS="\t"}
    NR==1 {print "chr", "pos", "rsid"}
    NR > 1 {
        rsid = $1 ":" $2 ":" $3 ":" $4
        print $1, $2, rsid
    }' some_gwas.tsv.tsv > variants.tsv
```

Six universes are built (one per gene set category) using `scripts/3_1_panukb_universes.sh`.

### 3.5 Bonferroni Threshold

P-value threshold for lead SNP selection:
```
0.05 / 6,685,228 = 7.479e-9
```
Where 6,685,228 is the number of variants in the reference panel (1000G EUR BIM).

---

## 4. Server Directory Structure

```
/media/DATA/gwasim/round2/
    bioGWAS/
        biogwas.py                          # Simulator source
        tests/
            data/                           # Reference data
                merged_1000genomes_eur.*    # 1000G EUR panel
                gencode.v37.annotation.gtf  # Gene annotations
                c2.cp.kegg.v2023.1.Hs.symbols.gmt  # C2 KEGG only (186 sets, for bioGWAS)
                c2.all.v2023.2.Hs.symbols.gmt      # C2 all curated (7233 sets, for Pan-UKB)
                c5.go.bp.v2024.1.Hs.symbols.gmt    # GO Biological Process (7608 sets)
                c5.go.cc.v2024.1.Hs.symbols.gmt    # GO Cellular Component (1026 sets)
                c5.go.mf.v2024.1.Hs.symbols.gmt    # GO Molecular Function (1820 sets)
                path_{small,medium,big,random}.txt  # Pathway gene lists
            3_pathways/
                in_data/                    # Continuous iter 0-29
                extra_in_data/              # Continuous iter 30-49 + sensitivity
                binary_in_data/             # Binary iter 0-49
    lsea_test/
        in_data/
            variants.tsv, anno.bed, uni.json, *.gmt
        lsea_results/                       # Original LSEA results (535 dirs)
        aggregated_data/                    # TPR/FPR CSVs for plotting
        validation_NEW_full/                # Re-run results (535 experiments)
    panukb/
        ukb_summstats/                      # Raw (*.tsv.tsv) + normalized (*.norm.tsv)
        download_panukb_summstats.py
        preproc_tsv.py
        run_preproc.sh
    panukb_lsea/
        in_data/
            variants.tsv, anno.bed, uni_*.json  # 6 universe files
        lsea_results/                       # Pan-UKB enrichment results (900 dirs)
        ukb_universe_bed/                   # BED files for GWAS-on-GWAS
        ukb_universe/                       # GWAS-on-GWAS universe
        aggregated_checks/                  # Summary TSVs per category
        validation_NEW_panukb/              # Re-run results
    LSEA/                                   # This repository
/media/DATA/bioinformatics/LSEA/
    tissues/
        GTEx8_formatted.gmt                 # GTEx tissue gene sets (45 non-empty tissues)
        blood_cell_markers.gmt              # Blood cell marker gene sets (12 sets)
```

---

## 5. Reproducing from Scratch

### bioGWAS (simulated data)
1. Install [bioGWAS](https://github.com/TohaRhymes/GWAS_simulator) Docker image
2. Prepare reference data (1000G EUR, GENCODE v37, GMT files)
3. Run simulation scripts (`1_iterate_pathways.py`, `extra_1_iterate_pathways.py`,
   `binary_1_iterate_pathways.py`, `k_1_iterate_pathways.py`)
4. Run `scripts/1_1_create_universe.sh` to build the universe
5. Run `scripts/1_2_lsea_continuous.sh`, `scripts/1_3_lsea_binary.sh`, `scripts/2_1_lsea_sensitivity.sh`

### Pan-UKB (real GWAS data)
1. Download Pan-UKB manifests and select phenotypes (`_select_panukb_h2.py`)
2. Download summary statistics (`download_panukb_summstats.py`)
3. Decompress and extract columns (`mv_to_dir_and_decompress_bgz.sh`)
4. Normalize for LSEA (`run_preproc.sh` using `preproc_tsv.py`)
5. Run `scripts/3_1_panukb_universes.sh` to build 6 universes
6. Run `scripts/3_2_panukb_enrichment.sh` (Experiment 1)
7. Run `scripts/4_1_panukb_gwas_on_gwas.sh` (Experiment 2)

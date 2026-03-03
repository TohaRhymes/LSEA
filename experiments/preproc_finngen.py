#!/usr/bin/env python3
"""
Preprocess FinnGen R9 summary statistics for LSEA/MAGMA/PASCAL analysis.

Steps:
1. Read gzipped FinnGen sumstats
2. Create BED for liftOver (GRCh38 → GRCh37)
3. Run liftOver
4. Merge lifted positions, create chr:pos:ref:alt variant IDs (GRCh37)
5. Filter to variants present in 1000G EUR .bim
6. Output:
   - O15_HYPTENSPREG.norm.tsv       (LSEA format: chr pos rsid pval)
   - O15_HYPTENSPREG_magma_pval.tsv (MAGMA format: rsid pval)
   - O15_HYPTENSPREG_pascal.txt     (PASCAL format: rsid pval, rs-prefixed only)
"""

import gzip
import subprocess
import sys
import os


def main(finngen_gz, chain_file, bim_file, liftover_bin, out_prefix):
    print(f"[INFO] Reading FinnGen sumstats: {finngen_gz}")

    # --- Step 1: Read FinnGen and create BED for liftOver ---
    variants = []  # (chrom, pos_hg38, ref, alt, pval, rsids)
    bed_path = out_prefix + "_liftover_input.bed"
    autosomal = set(str(c) for c in range(1, 23))

    need_bed = not os.path.exists(bed_path) or os.path.getsize(bed_path) == 0
    bed_out = open(bed_path, "w") if need_bed else None

    with gzip.open(finngen_gz, "rt") as f:
        header = f.readline().strip().lstrip("#").split("\t")
        idx = {col: i for i, col in enumerate(header)}

        required = ["chrom", "pos", "ref", "alt", "pval"]
        for col in required:
            if col not in idx:
                print(f"[ERROR] Column '{col}' not found in header: {header}")
                sys.exit(1)

        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[idx["chrom"]]
            pos = int(fields[idx["pos"]])
            ref = fields[idx["ref"]]
            alt = fields[idx["alt"]]
            pval = fields[idx["pval"]]
            rsids = fields[idx["rsids"]] if "rsids" in idx else ""

            # Skip non-autosomal and invalid
            if chrom not in autosomal:
                continue
            try:
                pval_f = float(pval)
                if pval_f <= 0 or pval_f > 1:
                    continue
            except ValueError:
                continue

            # unique ID for liftOver tracking
            var_id = f"{chrom}:{pos}:{ref}:{alt}"
            variants.append((chrom, pos, ref, alt, pval, rsids))

            # BED is 0-based: chr start end name
            if bed_out:
                bed_out.write(f"chr{chrom}\t{pos - 1}\t{pos}\t{var_id}\n")

    if bed_out:
        bed_out.close()

    print(f"[INFO] Read {len(variants)} autosomal variants from FinnGen")

    # --- Step 2: Run liftOver (skip if mapped file already exists) ---
    mapped_bed = out_prefix + "_mapped.bed"
    unmapped_bed = out_prefix + "_unmapped.bed"

    if os.path.exists(mapped_bed) and os.path.getsize(mapped_bed) > 0:
        print(f"[INFO] liftOver output already exists: {mapped_bed}, skipping liftOver")
    else:
        cmd = [liftover_bin, bed_path, chain_file, mapped_bed, unmapped_bed]
        print(f"[INFO] Running liftOver: {' '.join(cmd)}")
        ret = subprocess.call(cmd)
        if ret != 0:
            print(f"[ERROR] liftOver failed with exit code {ret}")
            sys.exit(1)

    # Count unmapped
    unmapped_count = 0
    if os.path.exists(unmapped_bed):
        with open(unmapped_bed) as f:
            for line in f:
                if not line.startswith("#"):
                    unmapped_count += 1
    print(f"[INFO] Unmapped variants: {unmapped_count}")

    # --- Step 3: Parse lifted positions ---
    # Map: original var_id -> (chrom_hg19, pos_hg19)
    lifted = {}
    with open(mapped_bed) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom_hg19 = fields[0].replace("chr", "")
            pos_hg19 = int(fields[2])  # end = 1-based pos (BED 0-based start, we want end)
            var_id = fields[3]
            lifted[var_id] = (chrom_hg19, pos_hg19)

    print(f"[INFO] Successfully lifted {len(lifted)} variants")

    # --- Step 4: Load 1000G EUR bim for filtering ---
    print(f"[INFO] Loading BIM file: {bim_file}")
    bim_rsids = set()
    with open(bim_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            bim_rsids.add(fields[1])  # col2 = variant ID (chr:pos:ref:alt)

    print(f"[INFO] BIM contains {len(bim_rsids)} variants")

    # --- Step 5: Merge and filter ---
    lsea_path = out_prefix + ".norm.tsv"
    magma_path = out_prefix + "_magma_pval.tsv"
    pascal_path = out_prefix + "_pascal.txt"

    n_matched = 0
    n_pascal = 0

    with open(lsea_path, "w") as f_lsea, \
         open(magma_path, "w") as f_magma, \
         open(pascal_path, "w") as f_pascal:

        f_lsea.write("chr\tpos\trsid\tpval\n")
        f_magma.write("rsid\tpval\n")

        for chrom, pos_hg38, ref, alt, pval, rsids in variants:
            var_id_hg38 = f"{chrom}:{pos_hg38}:{ref}:{alt}"

            if var_id_hg38 not in lifted:
                continue

            chrom_hg19, pos_hg19 = lifted[var_id_hg38]

            # Skip if chromosome changed during liftover
            if chrom_hg19 != chrom:
                continue

            # Create GRCh37 variant ID
            rsid_hg19 = f"{chrom}:{pos_hg19}:{ref}:{alt}"

            # Filter to 1000G EUR
            if rsid_hg19 not in bim_rsids:
                # Try flipped alleles (ref/alt swap)
                rsid_hg19_flip = f"{chrom}:{pos_hg19}:{alt}:{ref}"
                if rsid_hg19_flip in bim_rsids:
                    rsid_hg19 = rsid_hg19_flip
                else:
                    continue

            n_matched += 1

            # LSEA format
            f_lsea.write(f"{chrom}\t{pos_hg19}\t{rsid_hg19}\t{pval}\n")

            # MAGMA format (same rsid)
            f_magma.write(f"{rsid_hg19}\t{pval}\n")

            # PASCAL format: needs rs-prefixed IDs
            if rsids and rsids != "NA" and rsids != "":
                # FinnGen rsids can be comma-separated; take first rs-prefixed one
                for rs in rsids.split(","):
                    rs = rs.strip()
                    if rs.startswith("rs"):
                        f_pascal.write(f"{rs}\t{pval}\n")
                        n_pascal += 1
                        break

    print(f"[INFO] Matched to 1000G EUR: {n_matched} variants")
    print(f"[INFO] PASCAL (rs-prefixed): {n_pascal} variants")
    print(f"[INFO] Output files:")
    print(f"  LSEA:   {lsea_path}")
    print(f"  MAGMA:  {magma_path}")
    print(f"  PASCAL: {pascal_path}")

    # --- Cleanup intermediate files ---
    for f in [bed_path, mapped_bed, unmapped_bed]:
        if os.path.exists(f):
            os.remove(f)

    print("[INFO] Done!")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: preproc_finngen.py <finngen.gz> <chain_file> <bim_file> <liftover_bin> <out_prefix>")
        print("Example: preproc_finngen.py finngen_R9_O15_HYPTENSPREG.gz hg38ToHg19.over.chain merged_1000genomes_eur.bim ./liftOver O15_HYPTENSPREG")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

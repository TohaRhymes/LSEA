#!/usr/bin/env python3
"""
Preprocesses Pan-UKB GWAS summary statistics for LSEA.

Transformations:
  1. Converts neglog10_pval_EUR -> pval (= 10^(-neglog10_pval_EUR))
  2. Creates rsid column as chr:pos:ref:alt
  3. Filters to variants present in the PLINK BIM file

Usage:
    python3 preproc_tsv.py INPUT_FILE OUTPUT_FILE BIM_FILE

Example:
    python3 preproc_tsv.py \
        biomarkers-30600-both_sexes-irnt.tsv.tsv \
        biomarkers-30600-both_sexes-irnt.norm.tsv \
        /path/to/merged_1000genomes_eur.bim
"""

import pandas as pd
import sys


def process_data(input_file, output_file, rsid_file):
    bim = pd.read_csv(rsid_file, header=None, sep='\t')
    rsids = bim.iloc[:, 1]

    df = pd.read_csv(input_file, sep='\t')

    # Calculate the p-value from the negative log10 p-value
    df['pval'] = 10 ** (-df['neglog10_pval_EUR'])

    # Create the rsid column
    df['rsid'] = df['chr'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['ref'] + ':' + df['alt']

    df.drop(columns=['neglog10_pval_EUR'], inplace=True)

    # Filter to only variants present in the BIM file
    df[df.rsid.isin(rsids)].to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: ./preproc_tsv.py input_file output_file rsid_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    rsid_file = sys.argv[3]
    print("PREPROC Running:", input_file, output_file, rsid_file)
    process_data(input_file, output_file, rsid_file)

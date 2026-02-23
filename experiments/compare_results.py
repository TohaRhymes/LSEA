#!/usr/bin/env python3
"""
Compare annotation_stats TSV files between old and new LSEA results.

Usage:
    python3 experiments/compare_results.py OLD_DIR NEW_DIR

Example:
    python3 experiments/compare_results.py \
        /media/DATA/gwasim/round2/lsea_test/lsea_results \
        /media/DATA/gwasim/round2/lsea_test/validation_NEW_full
"""

import os
import sys
import csv
from pathlib import Path


def find_stats_files(base_dir):
    """Find all annotation_stats_*.tsv files under base_dir."""
    result = {}
    base = Path(base_dir)
    for stats_file in base.rglob("annotation_stats_*.tsv"):
        # Key: relative path from base_dir (e.g., "test10000_.../annotation_stats_uni.tsv")
        rel = stats_file.relative_to(base)
        result[str(rel)] = stats_file
    return result


def read_stats(path):
    """Read a stats TSV file and return rows as list of lists."""
    with open(path, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        return [row for row in reader]


def compare_stats(old_path, new_path):
    """Compare two stats files. Returns (identical, details)."""
    old_rows = read_stats(old_path)
    new_rows = read_stats(new_path)

    if old_rows == new_rows:
        return True, "identical"

    if len(old_rows) != len(new_rows):
        return False, f"row count differs: {len(old_rows)} vs {len(new_rows)}"

    diffs = []
    for i, (old_row, new_row) in enumerate(zip(old_rows, new_rows)):
        if old_row != new_row:
            diffs.append(f"  row {i}: OLD={old_row} NEW={new_row}")
    return False, "\n".join(diffs)


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    old_dir = sys.argv[1]
    new_dir = sys.argv[2]

    if not os.path.isdir(old_dir):
        print(f"ERROR: Old directory not found: {old_dir}")
        sys.exit(1)
    if not os.path.isdir(new_dir):
        print(f"ERROR: New directory not found: {new_dir}")
        sys.exit(1)

    old_files = find_stats_files(old_dir)
    new_files = find_stats_files(new_dir)

    all_keys = sorted(set(old_files.keys()) | set(new_files.keys()))

    identical = 0
    different = 0
    missing_old = 0
    missing_new = 0

    for key in all_keys:
        if key not in old_files:
            print(f"NEW_ONLY: {key}")
            missing_old += 1
        elif key not in new_files:
            print(f"OLD_ONLY: {key}")
            missing_new += 1
        else:
            ok, details = compare_stats(old_files[key], new_files[key])
            if ok:
                identical += 1
            else:
                print(f"DIFF: {key}")
                print(details)
                different += 1

    print()
    print("=" * 50)
    print(f"Identical:  {identical}")
    print(f"Different:  {different}")
    print(f"Old only:   {missing_new}")
    print(f"New only:   {missing_old}")
    print(f"Total keys: {len(all_keys)}")

    if different == 0 and missing_new == 0:
        print("\nAll matched — results are reproducible.")
    else:
        print("\nDifferences found — investigate before proceeding.")
        sys.exit(1)


if __name__ == "__main__":
    main()

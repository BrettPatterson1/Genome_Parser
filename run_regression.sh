#!/usr/bin/env bash

activate () {
  . ././venv/Scripts//activate
}

activate

python find_peptides.py GCF_000007465.2_ASM746v2_genomic.fna -a regression_anchor_locations.txt -r 500 -e test_result

diff regression_correct.csv test_result.csv

rm test_result.csv

echo "If no diff above, Passed"

read -n 1 -s -r -p "Press any key to continue"
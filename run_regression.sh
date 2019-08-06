#!/usr/bin/env bash

activate () {
  . ././venv/Scripts//activate
}

activate

python find_peptides.py GCF_000007465.2_ASM746v2_genomic.fna -e test_result

diff regression_correct.csv test_result.csv

rm test_result.csv

read -n 1 -s -r -p "Press any key to continue"
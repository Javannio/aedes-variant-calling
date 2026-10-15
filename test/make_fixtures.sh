#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data ref
for s in test_1 test_2; do
  for r in 1 2; do
    printf '@%s_r%s\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$s" "$r" \
      | gzip > "data/${s}_R${r}.fastq.gz"
  done
done
printf '>tiny_contig\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' \
  > ref/tiny.fa
ls -la data ref

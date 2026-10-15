#!/bin/bash
# Generate tiny paired-end test fixtures with planted SNPs so the CI run
# exercises real alignment + variant calling end-to-end.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data ref

python3 - <<'PY'
import gzip, random

random.seed(42)
ref = "".join(random.choice("ACGT") for _ in range(300))
with open("ref/tiny.fa", "w") as f:
    f.write(">tiny_contig\n")
    for i in range(0, len(ref), 60):
        f.write(ref[i:i+60] + "\n")

comp = str.maketrans("ACGT", "TGCA")

def revcomp(s):
    return s.translate(comp)[::-1]

READ_LEN = 50
INSERT = 120

def make_sample(name, n_pairs, snp_pos, snp_base):
    r1_path = f"data/{name}_R1.fastq.gz"
    r2_path = f"data/{name}_R2.fastq.gz"
    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for i in range(n_pairs):
            start = random.randint(0, len(ref) - INSERT - 1)
            frag = list(ref[start:start + INSERT])
            # plant heterozygous SNP in a subset of fragments
            if snp_pos is not None and start <= snp_pos < start + INSERT and i % 2 == 0:
                frag[snp_pos - start] = snp_base
            frag = "".join(frag)
            read1 = frag[:READ_LEN]
            read2 = revcomp(frag[-READ_LEN:])
            qual = "I" * READ_LEN
            r1.write(f"@{name}_{i}/1\n{read1}\n+\n{qual}\n")
            r2.write(f"@{name}_{i}/2\n{read2}\n+\n{qual}\n")

make_sample("test_1", 40, 150, "G" if ref[150] != "G" else "A")
make_sample("test_2", 40, 150, "G" if ref[150] != "G" else "A")
print(f"ref length: {len(ref)}; pos150={ref[150]}")
PY

ls -la data ref

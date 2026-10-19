#!/usr/bin/env python3
"""Generate inputs for the pareto optimization tests.

Usage: gen_inputs.py <genome fasta> <output directory>
small.fa, mid.fa -- prefixes of the genome (about 40k and 400k k-mers)
random.fa        -- sparse random 31-mers (no overlaps, many runs, small index type)
prefix_a.fa      -- k-mers sharing a prefix (empty buckets of the search index)
mixed.fa(.gz)    -- multi-line FASTA with N's and lowercase letters, also gzipped
"""
import gzip
import os
import random
import sys


def write_fasta(path, sequences, width=None):
    with open(path, "w") as f:
        for i, sequence in enumerate(sequences):
            f.write(f">{i}\n")
            if width:
                f.write("\n".join(sequence[j:j + width] for j in range(0, len(sequence), width)) + "\n")
            else:
                f.write(sequence + "\n")


def main():
    genome_path, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    genome = "".join(line.strip() for line in open(genome_path) if not line.startswith(">"))
    write_fasta(f"{out}/small.fa", [genome[:40000]])
    write_fasta(f"{out}/mid.fa", [genome[:400000]])

    rng = random.Random(1)
    write_fasta(f"{out}/random.fa", ["".join(rng.choice("ACGT") for _ in range(31)) for _ in range(3000)])
    write_fasta(f"{out}/prefix_a.fa", ["A" + "".join(rng.choice("AC") for _ in range(30)) for _ in range(300)])

    mixed = "".join(rng.choice("ACGTacgtN") if rng.random() < 0.05 else rng.choice("ACGT") for _ in range(20000))
    write_fasta(f"{out}/mixed.fa", [mixed, mixed[:5000].lower()], width=60)
    with open(f"{out}/mixed.fa", "rb") as source, gzip.open(f"{out}/mixed.fa.gz", "wb") as target:
        target.write(source.read())


if __name__ == "__main__":
    main()

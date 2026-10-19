#!/usr/bin/env python3
"""Verify that a masked superstring represents exactly the k-mers of a FASTA file.

Usage: verify_ms.py <fasta> <masked superstring> <k> [-u]
Prints the result with |S|, the number of runs of ones and of zeros in the mask, exits with 1 if the sets differ.
Without -u, k-mers and their reverse complements are treated as identical.
"""
import gzip
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical(kmer, bidirectional):
    return min(kmer, kmer.translate(COMPLEMENT)[::-1]) if bidirectional else kmer


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    sequences, current = [], []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                current = []
            else:
                current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def main():
    fasta, ms_path, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
    bidirectional = not (len(sys.argv) > 4 and sys.argv[4] == "-u")

    expected = set()
    for sequence in read_fasta(fasta):
        sequence = sequence.upper()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if set(kmer) <= set("ACGT"):
                expected.add(canonical(kmer, bidirectional))

    ms = "".join(read_fasta(ms_path))
    superstring = ms.upper()
    mask = [c.isupper() for c in ms]
    represented = {canonical(superstring[i:i + k], bidirectional)
                   for i in range(len(superstring) - k + 1) if mask[i]}

    runs = sum(1 for i, m in enumerate(mask) if m and (i == 0 or not mask[i - 1]))
    zeros = mask.count(False)
    ok = represented == expected
    print(f"ok={ok} |S|={len(superstring)} runs={runs} zeros={zeros} "
          f"missing={len(expected - represented)} extra={len(represented - expected)} n={len(expected)}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

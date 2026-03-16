#!/usr/bin/env python3

from sys import stdin, stderr, argv

def kmer_counter_ms(k: int):
    """
    Simply counts all distinct kmers in the input masked superstring.
    Very inefficiently. I have somehow lost the script from Pavel. XD
    """

    COMPLEMENTS = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }

    def get_canonical_kmer(kmer: str) -> str:
        complement = "".join([COMPLEMENTS[c] for c in reversed(kmer)])
        return min(kmer, complement)

    all = set()
    for line in stdin.readlines():
        line = line.strip()
        if line.startswith(">"): continue
        for i in range(len(line) - k + 1):
            if line[i] >= 'a': continue
            all.add(get_canonical_kmer(line[i:i+k].upper()))
    print(len(all))

if __name__ == "__main__":
    if len(argv) != 2:
        print("Usage: kmer_counter_ms.py <k> < input_file > <number_of_distinct_kmers>", file=stderr)
        exit(1)
    kmer_counter_ms(int(argv[1]))

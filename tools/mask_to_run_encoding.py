#!/usr/bin/env python3

from sys import argv

def matchtig_mask_to_sequence_lengths():
    """
    Prints the lengths of sequences corresponding to runs of '1's in the mask and K, separated by newlines.
    Works only with mask of simplitigs or matchtigs coverted to masked superstring --
    the assumption that every run of zeros has the same length k-1.
    We also need to print the K (length of runs of zeros, we print it as the last number), otherwise the information would be missing.
    Not very efficient implementation, but doesn't need to be at all.
    """
    last_char = False
    zeros_now = 0
    k = 1
    first_one = 0
    for i, c in enumerate(input()):
        if c == '0':
            if last_char: print(i - first_one)
            last_char = False
            zeros_now += 1
        elif c == '1':
            if not last_char:
                first_one = i
                if k == 1: k = zeros_now + 1
                else: assert k == zeros_now + 1, f"Matchtig mask to sequence length: Detected run of zeros with length {zeros_now} at position {i + 1}, expected {k - 1}."
                zeros_now = 0
            last_char = True
        else: raise ValueError(f"Matchtig mask to sequence length: Invalid character {c} at position {i}.")
    assert k > 0, f"Matchtig mask to sequence length: No run of zeros found in the mask, would mean K = 0."
    assert k == zeros_now + 1, f"Matchtig mask to sequence length: Detected run of zeros with length {zeros_now} at the end, expected {k - 1}."
    print(k)

def mask_to_sequence_ends(only_ones: bool = False):
    """
    Prints the positions of ends of runs of '1's and '0's in the mask, separated by newlines.
    If only_ones is set to True, only outputs positions of ends of runs of ones.
    Note that not using ends of runs of zeros means we lose the information about k, which cannot be simply appended for use in Elias-Fano encoding
    (it has to be accounted for separately; however, it usually only needs a byte).
    """
    last_char = True
    inp = input()
    for i, c in enumerate(inp):
        if (c == '1') != last_char:
            if not (only_ones and last_char): print(i)
            last_char = not last_char
    print(len(inp))

RPSS_METHODS = "unitigs|simplitigs|eulertigs|matchtigs|optimal_matchtigs"

if __name__ == "__main__":
    if len(argv) <= 1: matchtig_mask_to_sequence_lengths()
    elif len(argv) == 2:
        if argv[1] in RPSS_METHODS: mask_to_sequence_ends(True)
        else: mask_to_sequence_ends()
    else: print(f"Invalid number of arguments: {len(argv) - 1}. Expected 0 or 1.")

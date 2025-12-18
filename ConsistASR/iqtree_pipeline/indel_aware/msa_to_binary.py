#!/usr/bin/env python3
# ============================================================
# msa_to_binary.py
#
# Convert an MSA in FASTA format to a 0/1 (binary) representation
# for RAxML / RAxML-NG.
#
# Specification:
#   - gap ('-')      → 0
#   - any residue    → 1  (A–Z, X, etc. are all treated as "1")
#   - fully gapped columns are NOT removed
#   - output is a simple PHYLIP-like format
#
# Notes:
#   - Only the first token after ">" is used as the sequence name
#     (anything after the first whitespace is ignored).
#
# Usage:
#   python msa_to_binary.py input.fasta output_binary.phy
# ============================================================

import sys
from collections import OrderedDict


def read_fasta(filename):
    """
    Read a FASTA file into an OrderedDict: {name: sequence}.
    Only the first token on the header line is used as the name.
    """
    seqs = OrderedDict()
    name = None
    seq = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # flush previous sequence
                if name is not None:
                    seqs[name] = "".join(seq)
                # take first token as sequence name
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line)

    if name is not None:
        seqs[name] = "".join(seq)

    return seqs


def to_binary(msa_dict):
    """
    Convert an MSA dictionary {name: amino_acid_sequence}
    into a binary dictionary {name: 0/1 sequence}.
    """
    binary_dict = {}
    for name, seq in msa_dict.items():
        binary_seq = []
        for c in seq:
            if c == "-":
                binary_seq.append("0")   # gap = 0
            else:
                binary_seq.append("1")   # residue = 1
        binary_dict[name] = "".join(binary_seq)
    return binary_dict


def write_phylip(binary_dict, outfile):
    """
    Write a simple PHYLIP-like file:
        <nseq> <seqlen>
        <name> <binary_string>
    Names are left-padded up to 20 characters.
    """
    nseq = len(binary_dict)
    if nseq == 0:
        raise ValueError("No sequences to write (binary_dict is empty).")

    seqlen_set = {len(s) for s in binary_dict.values()}
    if len(seqlen_set) != 1:
        raise ValueError("Sequences have unequal lengths in binary_dict.")
    seqlen = seqlen_set.pop()

    with open(outfile, "w") as f:
        f.write(f"{nseq} {seqlen}\n")
        for name, seq in binary_dict.items():
            f.write(f"{name.ljust(20)} {seq}\n")


def main():
    if len(sys.argv) < 3:
        print("Usage: python msa_to_binary.py <input.fasta> <output.phy>", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    print(f"[INFO] Reading {infile} ...")
    msa = read_fasta(infile)
    nseq = len(msa)
    print(f"[INFO] {nseq} sequences loaded")

    if nseq == 0:
        print("[ERROR] No sequences were found in the input FASTA.", file=sys.stderr)
        sys.exit(1)

    # quick check for duplicate names
    if len(msa.keys()) != len(set(msa.keys())):
        print("[WARN] Duplicate sequence names detected in the FASTA headers.", file=sys.stderr)

    binary = to_binary(msa)
    length_set = {len(s) for s in binary.values()}
    if len(length_set) != 1:
        print("[WARN] Unequal sequence lengths detected in the alignment!", file=sys.stderr)

    write_phylip(binary, outfile)
    aln_len = list(length_set)[0]
    print(f"[DONE] Written to {outfile}")
    print(f"[INFO] Alignment length: {aln_len}")


if __name__ == "__main__":
    main()

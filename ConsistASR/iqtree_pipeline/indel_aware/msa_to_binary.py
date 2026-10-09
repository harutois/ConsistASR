#!/usr/bin/env python3
# ============================================================
# msa_to_binary.py
#
# Convert an amino-acid MSA in FASTA format to a 0/1 (binary)
# gap-state representation for RAxML / RAxML-NG.
#
# Specification:
#   - gap ('-')        → 0
#   - any non-gap char → 1
#   - fully gapped columns are NOT removed
#   - output is a simple PHYLIP-like format
#
# Notes:
#   - Only the first token after ">" is used as the sequence ID
#     (anything after the first whitespace is ignored).
#   - Sequence IDs must be unique after this first-token parsing.
#   - Input sequences must all have the same alignment length.
#   - Input should be a standard amino-acid MSA containing
#     amino-acid characters (including X if present) and '-' gaps.
#
# Usage:
#   python msa_to_binary.py input.fasta output_binary.phy
# ============================================================

import sys
from collections import OrderedDict


def read_fasta(filename):
    """
    Read a FASTA file into an OrderedDict: {name: sequence}.

    Only the first token on each FASTA header line is used as the
    sequence ID. Duplicate IDs are treated as errors.
    """
    seqs = OrderedDict()
    name = None
    seq = []

    def store_record(record_name, seq_chunks):
        """Store one FASTA record after checking for duplicate IDs."""
        if record_name in seqs:
            raise ValueError(
                f"Duplicate FASTA sequence ID detected: '{record_name}'. "
                "Sequence IDs must be unique."
            )
        seqs[record_name] = "".join(seq_chunks)

    with open(filename) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                # Store previous sequence
                if name is not None:
                    store_record(name, seq)

                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Empty FASTA header detected at line {lineno}."
                    )

                # Use only the first token as the sequence ID
                name = header.split()[0]
                seq = []

            else:
                if name is None:
                    raise ValueError(
                        f"Sequence data found before the first FASTA header "
                        f"at line {lineno}."
                    )
                seq.append(line)

    # Store final sequence
    if name is not None:
        store_record(name, seq)

    return seqs


def to_binary(msa_dict):
    """
    Convert an MSA dictionary {name: amino_acid_sequence}
    into a binary dictionary {name: 0/1 sequence}.

    '-' is encoded as 0; every other character is encoded as 1.
    """
    binary_dict = OrderedDict()

    for name, seq in msa_dict.items():
        binary_seq = []

        for c in seq:
            if c == "-":
                binary_seq.append("0")
            else:
                binary_seq.append("1")

        binary_dict[name] = "".join(binary_seq)

    return binary_dict


def write_phylip(binary_dict, outfile):
    """
    Write a simple PHYLIP-like file:

        <nseq> <seqlen>
        <name> <binary_string>

    Sequence names shorter than 20 characters are right-padded
    with spaces to a minimum width of 20 characters. Longer names
    are written without truncation.
    """
    nseq = len(binary_dict)

    if nseq == 0:
        raise ValueError("No sequences to write (binary_dict is empty).")

    seqlen_set = {len(seq) for seq in binary_dict.values()}

    if len(seqlen_set) != 1:
        raise ValueError(
            "Sequences have unequal lengths; input must be an aligned MSA."
        )

    seqlen = next(iter(seqlen_set))

    with open(outfile, "w") as f:
        f.write(f"{nseq} {seqlen}\n")

        for name, seq in binary_dict.items():
            f.write(f"{name.ljust(20)} {seq}\n")


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: python msa_to_binary.py <input.fasta> <output.phy>",
            file=sys.stderr,
        )
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    try:
        print(f"[INFO] Reading {infile} ...")
        msa = read_fasta(infile)

        nseq = len(msa)
        print(f"[INFO] {nseq} sequences loaded")

        if nseq == 0:
            raise ValueError("No sequences were found in the input FASTA.")

        # Check that the input really is an alignment
        length_set = {len(seq) for seq in msa.values()}

        if len(length_set) != 1:
            raise ValueError(
                "Unequal sequence lengths detected; "
                "input must be an aligned MSA."
            )

        aln_len = next(iter(length_set))

        binary = to_binary(msa)
        write_phylip(binary, outfile)

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"[DONE] Written to {outfile}")
    print(f"[INFO] Alignment length: {aln_len}")


if __name__ == "__main__":
    main()
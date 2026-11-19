#!/usr/bin/env python3
# ============================================================
# msa_to_binary.py
#
# Convert an amino-acid MSA in FASTA format to a 0/1 (binary)
# gap-state representation for RAxML / RAxML-NG.
#
# Specification:
#   - gap ('-')                         → 0
#   - standard amino acids             → 1
#   - ambiguity/unknown codes B/Z/J/X  → 1
#   - any other character              → ERROR
#   - fully gapped columns are NOT removed
#   - output is a simple PHYLIP-like format
#
# Notes:
#   - Only the first token after ">" is used as the sequence ID
#     (anything after the first whitespace is ignored).
#   - Sequence IDs must be unique after this first-token parsing.
#   - Input sequences must all have the same alignment length.
#   - Lowercase amino-acid letters are accepted and interpreted
#     equivalently to uppercase letters.
#   - Only '-' is interpreted as an alignment gap. Symbols such as
#     '.', '?', '~', '*', and '!' are rejected rather than silently
#     encoded as residue presence.
#
# Usage:
#   python msa_to_binary.py input.fasta output_binary.phy
# ============================================================

import sys
from collections import OrderedDict


# Standard amino acids plus common ambiguity/unknown residue codes.
# These symbols indicate that a residue is present at the alignment
# position even if its exact amino-acid identity is uncertain.
VALID_PRESENT = set("ACDEFGHIKLMNPQRSTVWYBZJX")
GAP_CHAR = "-"


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

    Encoding
    --------
    '-'                    -> 0
    standard AA/B/Z/J/X   -> 1
    any other character   -> error

    Lowercase amino-acid letters are accepted by converting them
    to uppercase for validation.
    """
    binary_dict = OrderedDict()

    for name, seq in msa_dict.items():
        binary_seq = []

        for position, raw_char in enumerate(seq, start=1):
            char = raw_char.upper()

            if char == GAP_CHAR:
                binary_seq.append("0")

            elif char in VALID_PRESENT:
                binary_seq.append("1")

            else:
                hint = ""

                if raw_char in {".", "?", "~", "*", "!"}:
                    hint = (
                        " Only '-' is accepted as an alignment gap; "
                        "missing/other special symbols are not assigned "
                        "a binary presence/absence state automatically."
                    )

                raise ValueError(
                    f"Unexpected character '{raw_char}' in sequence "
                    f"'{name}' at alignment position {position}."
                    f"{hint}"
                )

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
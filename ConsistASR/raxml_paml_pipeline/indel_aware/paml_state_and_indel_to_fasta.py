#!/usr/bin/env python3
# ============================================================
# paml_state_and_indel_to_fasta.py
#
# Combine:
#   - PAML rst ancestral amino-acid sequences, and
#   - RAxML binary gap-state profiles already mapped to
#     PAML node IDs
#
# to produce:
#   1) raw PAML ancestral FASTA at alignment length
#   2) raw gap-stripped PAML ancestral FASTA
#   3) indel-aware ancestral FASTA with gaps
#   4) indel-aware gap-stripped ancestral FASTA
# ============================================================

import argparse
import re
import sys


# ---------------------------
# Extract AA sequences from PAML rst
# ---------------------------

def parse_all_nodes_from_rst(rst_path):
    """
    Parse ancestral amino-acid sequences from the PAML rst
    'List of extant and reconstructed sequences' section.

    Returns
    -------
    dict[str, str]
        {node_id: amino_acid_sequence}
    """
    node_seqs = {}
    current_id = None
    in_sequence_section = False
    found_sequence_section = False

    node_re = re.compile(
        r"^\s*node\s*#\s*(\d+)\b(.*)$",
        flags=re.IGNORECASE,
    )

    def parse_sequence_chunk(text, context):
        """
        Remove whitespace from a sequence fragment and verify
        that it contains only amino-acid letters or '-'.
        """
        compact = re.sub(r"\s+", "", text)

        if not compact:
            return ""

        if not re.fullmatch(r"[A-Za-z-]+", compact):
            raise ValueError(
                f"Unexpected character(s) in PAML ancestral "
                f"sequence ({context}): '{text.strip()}'"
            )

        return compact.upper()

    with open(rst_path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            raw = line.rstrip("\r\n")
            stripped = raw.strip()
            lower = stripped.lower()

            # Enter reconstructed-sequence section
            if "list of extant and reconstructed sequences" in lower:
                in_sequence_section = True
                found_sequence_section = True
                current_id = None
                continue

            if not in_sequence_section:
                continue

            # End of sequence section
            if (
                lower.startswith("overall")
                or lower.startswith("prob distribution at node")
            ):
                break

            # Detect ancestral sequence line
            m = node_re.match(raw)

            if m:
                node_id = m.group(1)

                if node_id in node_seqs:
                    raise ValueError(
                        f"Duplicate ancestral node #{node_id} "
                        f"in PAML rst at line {lineno}."
                    )

                current_id = node_id
                node_seqs[current_id] = []

                tail = m.group(2)

                chunk = parse_sequence_chunk(
                    tail,
                    f"node {node_id}, line {lineno}",
                )

                if chunk:
                    node_seqs[current_id].append(chunk)

                continue

            # Continuation line for a reconstructed sequence
            if current_id is not None:
                if not stripped:
                    current_id = None
                    continue

                chunk = parse_sequence_chunk(
                    raw,
                    f"node {current_id}, line {lineno}",
                )

                if chunk:
                    node_seqs[current_id].append(chunk)

    if not found_sequence_section:
        raise ValueError(
            "'List of extant and reconstructed sequences' "
            f"section not found in {rst_path}"
        )

    seq_dict = {
        node_id: "".join(chunks)
        for node_id, chunks in node_seqs.items()
    }

    # Reject empty reconstructed sequences
    empty_nodes = [
        node_id
        for node_id, seq in seq_dict.items()
        if not seq
    ]

    if empty_nodes:
        raise ValueError(
            "Empty ancestral sequence(s) found for PAML node(s): "
            + ", ".join(sorted(empty_nodes, key=int))
        )

    if not seq_dict:
        raise ValueError(
            f"No ancestral node sequences found in {rst_path}"
        )

    # All ancestral sequences should represent the same alignment.
    lengths = {
        len(seq)
        for seq in seq_dict.values()
    }

    if len(lengths) != 1:
        raise ValueError(
            "PAML ancestral sequences have inconsistent lengths: "
            + ", ".join(str(x) for x in sorted(lengths))
        )

    alignment_length = next(iter(lengths))

    print(
        f"[INFO] PAML ancestral alignment length: "
        f"{alignment_length}"
    )

    return seq_dict


# ---------------------------
# Parse binary gap-state file
# ---------------------------

def parse_indel_file(indel_file):
    """
    Parse RAxML binary gap-state profiles already mapped
    to PAML node IDs.

    Expected format:
        <node_name>  1010010101...

    Returns
    -------
    dict[str, list[str]]
        {node_name: list_of_bits}
    """
    indels = {}

    with open(indel_file, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            stripped = line.strip()

            if not stripped:
                continue

            parts = stripped.split()

            if len(parts) < 2:
                raise ValueError(
                    f"Malformed binary gap-state line "
                    f"at line {lineno}: '{stripped}'"
                )

            name = parts[0]
            bits = parts[1]

            if name in indels:
                raise ValueError(
                    f"Duplicate node name in binary "
                    f"gap-state file: '{name}'"
                )

            if not re.fullmatch(r"[01]+", bits):
                raise ValueError(
                    f"Invalid binary gap-state pattern "
                    f"for node {name} at line {lineno}. "
                    "Only '0' and '1' are allowed."
                )

            indels[name] = list(bits)

    if not indels:
        raise ValueError(
            f"No binary gap-state profiles found in {indel_file}"
        )

    return indels


# ---------------------------
# Merge AA states with binary gap states
# ---------------------------

def merge_state_and_indel(state_dict, indel_dict):
    """
    Merge PAML ancestral amino-acid sequences with binary
    gap-state profiles.

    Rules
    -----
    bit = "0"
        Inferred gap/absence state -> write '-'

    bit = "1"
        Inferred residue/presence state -> retain the PAML
        amino-acid state.
    """
    merged = {}
    missing_indel = []

    for node, seq in state_dict.items():

        if node not in indel_dict:
            missing_indel.append(node)
            continue

        bits = indel_dict[node]

        if len(bits) != len(seq):
            raise ValueError(
                f"Length mismatch for node {node}: "
                f"state={len(seq)}, indel={len(bits)}"
            )

        new_seq = []

        for aa, bit in zip(seq, bits):

            if bit == "0":
                new_seq.append("-")

            elif bit == "1":
                new_seq.append(aa)

            else:
                # Defensive check; parse_indel_file() should
                # already prevent this.
                raise ValueError(
                    f"Unexpected binary state '{bit}' "
                    f"for node {node}."
                )

        merged[node] = "".join(new_seq)

    # Every PAML ancestral node must have a mapped binary profile.
    if missing_indel:
        raise ValueError(
            f"{len(missing_indel)} PAML ancestral node(s) "
            "are missing from the binary gap-state file: "
            + ", ".join(
                sorted(missing_indel, key=int)
            )
        )

    # Extra RAxML nodes can arise from internal-node resolution
    # differences. They are not PAML ancestral nodes and are ignored.
    extra_indel_nodes = (
        set(indel_dict)
        - set(state_dict)
    )

    if extra_indel_nodes:
        print(
            f"[WARN] {len(extra_indel_nodes)} node(s) found "
            "in binary gap-state file but not in PAML rst "
            "(ignored): "
            + ", ".join(
                sorted(extra_indel_nodes)
            ),
            file=sys.stderr,
        )

    return merged


# ---------------------------
# FASTA writer
# ---------------------------

def write_fasta(seqs, outpath, remove_gaps=False):
    """
    Write sequences to FASTA.

    If remove_gaps=True, '-' characters are removed before
    writing.
    """
    with open(
        outpath,
        "w",
        encoding="utf-8",
    ) as f:

        for name, seq in seqs.items():

            out_seq = (
                seq.replace("-", "")
                if remove_gaps
                else seq
            )

            f.write(
                f">{name}\n{out_seq}\n"
            )


# ---------------------------
# Main
# ---------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Combine PAML rst ancestral amino-acid sequences "
            "with RAxML binary gap-state profiles to produce "
            "raw and indel-aware ancestral FASTA files."
        )
    )

    parser.add_argument(
        "--rst",
        required=True,
        help="PAML rst file containing ancestral sequences",
    )

    parser.add_argument(
        "--indel",
        required=True,
        help=(
            "RAxML binary gap-state file with node labels "
            "already mapped to PAML node IDs"
        ),
    )

    parser.add_argument(
        "--out_raw_withgap",
        default=None,
        help="Output raw PAML ASR FASTA at alignment length",
    )

    parser.add_argument(
        "--out_raw_nogap",
        default=None,
        help="Output raw PAML ASR FASTA with '-' removed",
    )

    parser.add_argument(
        "--out_withgap",
        default="ASR_PAML_indel_withgap.fasta",
        help="Output indel-aware FASTA with gaps",
    )

    parser.add_argument(
        "--out_nogap",
        default="ASR_PAML_indel_nogap.fasta",
        help="Output gap-stripped indel-aware FASTA",
    )

    args = parser.parse_args()

    try:

        print(
            f"[INFO] Reading PAML rst: {args.rst}"
        )

        state_dict = parse_all_nodes_from_rst(
            args.rst
        )

        print(
            f"[INFO] Loaded {len(state_dict)} "
            "ancestral nodes from rst"
        )

        print(
            f"[INFO] Reading binary gap-state file: "
            f"{args.indel}"
        )

        indel_dict = parse_indel_file(
            args.indel
        )

        print(
            f"[INFO] Loaded {len(indel_dict)} "
            "binary gap-state profiles"
        )

        print(
            "[INFO] Merging PAML AA states and "
            "binary gap-state information..."
        )

        merged = merge_state_and_indel(
            state_dict,
            indel_dict,
        )

        print(
            f"[INFO] Merged {len(merged)} sequences"
        )

        # Defensive final consistency check
        if len(merged) != len(state_dict):
            raise ValueError(
                f"Merged sequence count ({len(merged)}) "
                "does not match the number of PAML "
                f"ancestral sequences ({len(state_dict)})."
            )

        if args.out_raw_withgap:
            write_fasta(
                state_dict,
                args.out_raw_withgap,
                remove_gaps=False,
            )

        if args.out_raw_nogap:
            write_fasta(
                state_dict,
                args.out_raw_nogap,
                remove_gaps=True,
            )

        write_fasta(
            merged,
            args.out_withgap,
            remove_gaps=False,
        )

        write_fasta(
            merged,
            args.out_nogap,
            remove_gaps=True,
        )

    except (OSError, ValueError) as exc:

        print(
            f"[ERROR] {exc}",
            file=sys.stderr,
        )

        sys.exit(1)

    print("[DONE] Output written:")

    if args.out_raw_withgap:
        print(
            f"           raw with gaps     : "
            f"{args.out_raw_withgap}"
        )

    if args.out_raw_nogap:
        print(
            f"           raw no gaps       : "
            f"{args.out_raw_nogap}"
        )

    print(
        f"           indel with gaps   : "
        f"{args.out_withgap}"
    )

    print(
        f"           indel no gaps     : "
        f"{args.out_nogap}"
    )


if __name__ == "__main__":
    main()
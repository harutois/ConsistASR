#!/usr/bin/env python3
# ============================================================
# state_and_indel_to_fasta.py
#
# Combine IQ-TREE .state output and RAxML-based binary gap-state
# results into indel-aware ancestral FASTA sequences.
#
# Outputs:
#   1) Raw IQ-TREE ancestral FASTA at alignment length
#   2) Raw gap-stripped IQ-TREE ancestral FASTA
#   3) Indel-aware alignment-like FASTA with gaps
#   4) Indel-aware gap-stripped FASTA
# ============================================================

import argparse
import sys
from collections import defaultdict


# ---------------------------
# Extract AA sequences for all nodes from IQ-TREE .state
# ---------------------------

def parse_state_file(state_file):
    """
    Parse an IQ-TREE .state file and return sequences per node.

    Expected columns:
        Node  Site  State  p_A  p_R ...

    Returns
    -------
    dict[str, str]
        {node_name: amino_acid_sequence}

    Notes
    -----
    Site numbers are checked for each node and must form a
    contiguous 1-based series (1, 2, 3, ...).
    """
    states = defaultdict(list)

    with open(state_file) as f:
        for lineno, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split()

            if len(parts) < 3:
                continue

            # Skip the IQ-TREE .state header:
            # Node  Site  State  p_A ...
            if (
                parts[0] == "Node"
                and parts[1] == "Site"
                and parts[2] == "State"
            ):
                continue

            node = parts[0]

            try:
                site = int(parts[1])
            except ValueError:
                raise ValueError(
                    f"Invalid site number in .state file at line {lineno}: "
                    f"'{parts[1]}'"
                )

            aa = parts[2]
            states[node].append((site, aa))

    if not states:
        raise ValueError(
            "No ancestral states were found in the IQ-TREE .state file."
        )

    seqs = {}

    for node, site_states in states.items():
        site_states.sort(key=lambda x: x[0])

        sites = [site for site, _ in site_states]
        expected_sites = list(range(1, len(site_states) + 1))

        if sites != expected_sites:
            raise ValueError(
                f"Non-contiguous or duplicated site numbers detected "
                f"for node {node}. Expected sites 1-{len(site_states)}, "
                f"but found an inconsistent site series."
            )

        seqs[node] = "".join(aa for _, aa in site_states)

    return seqs


# ---------------------------
# Parse binary gap-state file
# ---------------------------

def parse_indel_file(indel_file):
    """
    Parse the RAxML binary gap-state output after node-name mapping.

    Expected format:
        Node123  1010100111...

    Returns
    -------
    dict[str, list[str]]
        {node_name: list_of_bits}

    Notes
    -----
    Only '0' and '1' are accepted in the binary state string.
    Duplicate node names are treated as errors.
    """
    indels = {}

    with open(indel_file) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 2:
                raise ValueError(
                    f"Malformed indel-state line at line {lineno}: '{line}'"
                )

            name = parts[0]
            bits = parts[1].strip()

            if name in indels:
                raise ValueError(
                    f"Duplicate node name in indel file: '{name}'"
                )

            invalid_chars = set(bits) - {"0", "1"}
            if invalid_chars:
                invalid_display = "".join(sorted(invalid_chars))
                raise ValueError(
                    f"Invalid binary state for node {name}: "
                    f"unexpected character(s) '{invalid_display}'. "
                    "Only '0' and '1' are allowed."
                )

            indels[name] = list(bits)

    if not indels:
        raise ValueError(
            "No binary gap-state profiles were found in the indel file."
        )

    return indels


# ---------------------------
# Merge ASR AA states with binary gap-state patterns
# ---------------------------

def merge_state_and_indel(state_dict, indel_dict):
    """
    Merge amino-acid ASR states with binary gap-state information.

    Rules
    -----
    bit = "0"
        Inferred gap/absence state -> write '-'

    bit = "1"
        Inferred residue/presence state -> retain the amino-acid
        state reconstructed by IQ-TREE.

    Parameters
    ----------
    state_dict : dict[str, str]
        {node_name: amino_acid_sequence}

    indel_dict : dict[str, list[str]]
        {node_name: list_of_bits}

    Returns
    -------
    dict[str, str]
        {node_name: merged_sequence_with_gaps}
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
                # This should already have been caught by parse_indel_file(),
                # but retain the check here as a safeguard.
                raise ValueError(
                    f"Unexpected binary state '{bit}' for node {node}."
                )

        merged[node] = "".join(new_seq)

    # Missing binary profile for an IQ-TREE ancestral node is a hard error.
    # Otherwise, the final FASTA would silently lose that ancestral node.
    if missing_indel:
        raise ValueError(
            f"{len(missing_indel)} node(s) found in .state but missing "
            "from the indel file: "
            + ", ".join(sorted(missing_indel))
        )

    # Extra RAxML binary nodes can occur because of differences in
    # internal-node resolution near a polytomy or root. These are ignored.
    extra_indel_nodes = set(indel_dict) - set(state_dict)

    if extra_indel_nodes:
        print(
            f"[WARN] {len(extra_indel_nodes)} node(s) found in indel file "
            "but not in .state (ignored): "
            + ", ".join(sorted(extra_indel_nodes)),
            file=sys.stderr,
        )

    return merged


# ---------------------------
# FASTA writer
# ---------------------------

def write_fasta(seqs, outpath, remove_gaps=False):
    """
    Write sequences to FASTA.

    Parameters
    ----------
    seqs : dict[str, str]
        {name: sequence}

    outpath : str
        Output FASTA path.

    remove_gaps : bool, optional
        If True, remove all '-' characters before writing.
    """
    with open(outpath, "w") as f:
        for name, seq in seqs.items():
            out_seq = seq.replace("-", "") if remove_gaps else seq
            f.write(f">{name}\n{out_seq}\n")


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine IQ-TREE .state output and RAxML binary gap-state "
            "results into raw and indel-aware ancestral FASTA files."
        )
    )

    parser.add_argument(
        "--state",
        required=True,
        help="IQ-TREE .state file",
    )

    parser.add_argument(
        "--indel",
        required=True,
        help=(
            "RAxML binary gap-state file already mapped to "
            "IQ-TREE node names"
        ),
    )

    parser.add_argument(
        "--out_raw_withgap",
        default=None,
        help="Output raw IQ-TREE ASR FASTA at alignment length",
    )

    parser.add_argument(
        "--out_raw_nogap",
        default=None,
        help="Output raw IQ-TREE ASR FASTA with '-' removed",
    )

    parser.add_argument(
        "--out_withgap",
        default="ASR_final_withgap.fasta",
        help="Output indel-aware FASTA with gaps",
    )

    parser.add_argument(
        "--out_nogap",
        default="ASR_final_nogap.fasta",
        help="Output gap-stripped indel-aware FASTA",
    )

    args = parser.parse_args()

    try:
        print(f"[INFO] Reading state file: {args.state}")
        state_dict = parse_state_file(args.state)
        print(
            f"[INFO] Loaded {len(state_dict)} sequences from state file"
        )

        print(f"[INFO] Reading binary gap-state file: {args.indel}")
        indel_dict = parse_indel_file(args.indel)
        print(
            f"[INFO] Loaded {len(indel_dict)} binary gap-state profiles"
        )

        print("[INFO] Merging state and binary gap-state information...")
        merged = merge_state_and_indel(
            state_dict,
            indel_dict,
        )

        print(f"[INFO] Merged {len(merged)} sequences")

        # Defensive consistency check:
        # every IQ-TREE ancestral node should appear in the merged output.
        if len(merged) != len(state_dict):
            raise ValueError(
                f"Merged sequence count ({len(merged)}) does not match "
                f"the number of IQ-TREE ancestral sequences "
                f"({len(state_dict)})."
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
        print(f"[ERROR] {exc}", file=sys.stderr)
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
        f"           indel with gaps   : {args.out_withgap}"
    )
    print(
        f"           indel no gaps     : {args.out_nogap}"
    )


if __name__ == "__main__":
    main()
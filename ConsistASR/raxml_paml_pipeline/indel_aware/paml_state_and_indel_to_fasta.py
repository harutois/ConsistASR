#!/usr/bin/env python3
# ============================================================
# paml_state_and_indel_to_fasta.py
#
# Merge:
#   - a PAML .rst file (ASR results), and
#   - a RAxML indel 0/1 file (node labels already mapped to
#     the PAML node IDs)
# to produce:
#
#   - gap-containing FASTA with indels applied, and
#   - gap-stripped FASTA (e.g., for AlphaFold input)
#
# for all ancestral nodes.
#
# The interface and behavior are analogous to the IQ-TREE version
# (state_and_indel_to_fasta.py), while the .rst parsing logic is
# generalized from extract_multi_node_indel.py.
# ============================================================

import argparse
import re
from collections import defaultdict


# ---------------------------
# Extract AA sequences for all nodes from PAML .rst
# ---------------------------

def parse_all_nodes_from_rst(rst_path):
    """
    Parse all "node # <id>" blocks from a PAML .rst file and return:

        { '229': 'AA....', '230': 'AA....', ... }

    This is a generalization of parse_single_node_from_rst()
    from extract_multi_node_indel.py to handle all nodes at once.
    """
    node_seqs = defaultdict(list)
    current_id = None
    grabbing = False

    with open(rst_path, encoding="utf-8") as f:
        for line in f:
            raw = line.rstrip("\n")

            # Detect header line: "node # <id>"
            m = re.match(r"^\s*node\s+#\s*(\d+)\b", raw, flags=re.IGNORECASE)
            if m:
                node_id = m.group(1)  # e.g., '229'
                current_id = node_id
                grabbing = True

                # Some AA characters may already appear on the header line
                tail = raw.split("#", 1)[-1]
                pieces = re.findall(r"[A-Z]+", tail)
                if pieces:
                    node_seqs[current_id].append("".join(pieces))
                continue

            if grabbing:
                s = raw.strip()
                # Block termination: next node or TREE/Ancestral/Probab/empty line
                if (
                    s.startswith("node #")
                    or s.startswith("TREE")
                    or s.startswith("Ancestral")
                    or s.startswith("Probab")
                    or s == ""
                ):
                    current_id = None
                    grabbing = False
                    continue

                letters = re.findall(r"[A-Z]+", s)
                if letters and current_id is not None:
                    node_seqs[current_id].append("".join(letters))

    # Join collected chunks per node
    seq_dict = {nid: "".join(chunks) for nid, chunks in node_seqs.items()}

    if not seq_dict:
        raise ValueError(f"No node sequences found in {rst_path}")

    # Optional: length consistency check
    lengths = {len(s) for s in seq_dict.values()}
    if len(lengths) != 1:
        print(f"[WARN] PAML .rst sequences have multiple lengths: {lengths}")

    return seq_dict


# ---------------------------
# Parse indel (0/1 pattern) file
# ---------------------------

def parse_indel_file(indel_file):
    """
    Parse a RAxML-derived indel file.

    Expected format (one node per line):
        <node_name>  1010010101...

    Example:
        229  101010010101...

    Returns
    -------
    dict[str, list[str]]
        { node_name : ['1','0','1', ...] }
    """
    indels = {}
    with open(indel_file, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, bits = parts[0], parts[1].strip()
            if not re.fullmatch(r"[01]+", bits):
                # Ignore non-pattern lines (e.g., headers)
                continue
            indels[name] = list(bits)

    if not indels:
        raise ValueError(f"No valid 0/1 patterns found in {indel_file}")

    return indels


# ---------------------------
# Merge ASR AA states with indel patterns
# ---------------------------

def merge_state_and_indel(state_dict, indel_dict):
    """
    Merge AA sequences from PAML .rst with RAxML indel (0/1) patterns.

    Rules:
      bit = '0' → '-' (deletion) at that position
      bit = '1' → keep the AA as is (if AA is already '-', keep '-')

    Parameters
    ----------
    state_dict : dict[str, str]
        { node_name : AA_sequence }
        e.g., { '229': 'MKT...', ... }
    indel_dict : dict[str, list[str]]
        { node_name : list_of_bits }
        e.g., { '229': ['1','0',...], ... }

    Returns
    -------
    dict[str, str]
        { node_name : merged_sequence_with_indels }
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
                f"[ERROR] Length mismatch for node {node}: "
                f"state={len(seq)} indel={len(bits)}"
            )

        new_seq = []
        for aa, bit in zip(seq, bits):
            if bit == "0":
                new_seq.append("-")
            else:
                # bit == "1" → use AA as is (preserve '-' if present)
                new_seq.append(aa if aa != "-" else "-")
        merged[node] = "".join(new_seq)

    if missing_indel:
        print(
            f"[WARN] {len(missing_indel)} nodes found in .rst "
            f"but missing from indel file (skipped). "
            f"Examples: {missing_indel[:5]}"
        )

    extra_indel_nodes = set(indel_dict.keys()) - set(state_dict.keys())
    if extra_indel_nodes:
        print(
            f"[WARN] {len(extra_indel_nodes)} nodes found in indel file "
            f"but not in .rst (ignored). "
            f"Examples: {sorted(list(extra_indel_nodes))[:5]}"
        )

    return merged


# ---------------------------
# FASTA writer
# ---------------------------

def write_fasta(seqs, outpath, remove_gaps=False):
    """
    Write a dict of sequences {name: seq} to a FASTA file.

    If remove_gaps=True, all '-' characters are stripped from each sequence
    before writing (matching the IQ-TREE version behavior).
    """
    with open(outpath, "w", encoding="utf-8") as f:
        for name, seq in seqs.items():
            out_seq = seq.replace("-", "") if remove_gaps else seq
            f.write(f">{name}\n")
            # For now, write the sequence as a single line (no wrapping).
            f.write(out_seq + "\n")


# ---------------------------
# main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine PAML .rst ASR output and RAxML-based indel (0/1) results "
            "into gap-aware ASR FASTA (with and without gaps) for all nodes."
        )
    )
    parser.add_argument("--rst", required=True,
                        help="PAML .rst file (ASR result)")
    parser.add_argument("--indel", required=True,
                        help="RAxML indel file (0/1 patterns, node labels already mapped to PAML node IDs)")
    parser.add_argument("--out_withgap", default="ASR_PAML_indel_withgap.fasta",
                        help="Output FASTA (alignment-like, with gaps)")
    parser.add_argument("--out_nogap", default="ASR_PAML_indel_nogap.fasta",
                        help="Output FASTA (gap-stripped, e.g., for AlphaFold input)")
    args = parser.parse_args()

    print(f"[INFO] Reading PAML .rst: {args.rst}")
    state_dict = parse_all_nodes_from_rst(args.rst)
    print(f"[INFO] Loaded {len(state_dict)} nodes from .rst")

    print(f"[INFO] Reading indel file: {args.indel}")
    indel_dict = parse_indel_file(args.indel)
    print(f"[INFO] Loaded {len(indel_dict)} indel profiles")

    print("[INFO] Merging AA states and indel patterns...")
    merged = merge_state_and_indel(state_dict, indel_dict)
    print(f"[INFO] Merged {len(merged)} sequences")

    write_fasta(merged, args.out_withgap, remove_gaps=False)
    write_fasta(merged, args.out_nogap, remove_gaps=True)

    print("[DONE] Output written:")
    print(f"           with gaps : {args.out_withgap}")
    print(f"           no gaps   : {args.out_nogap}")


if __name__ == "__main__":
    main()

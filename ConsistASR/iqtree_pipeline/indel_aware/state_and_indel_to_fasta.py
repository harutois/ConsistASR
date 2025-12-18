#!/usr/bin/env python3
# ============================================================
# state_and_indel_to_fasta.py
#
# Combine IQ-TREE .state output and RAxML-based indel results
# (0/1 matrix) into indel-aware ancestral FASTA sequences.
#
# Outputs:
#   1) Alignment-like FASTA (with gaps preserved)
#   2) Gap-stripped FASTA (for AlphaFold input, etc.)
# ============================================================

import argparse
from collections import defaultdict


def parse_state_file(state_file):
    """
    Parse an IQ-TREE .state file and return sequences per node.

    Assumes each non-comment line has at least:
        <node> <site_index> <AA> ...
    and concatenates column 3 (AA) along the alignment.

    Returns
    -------
    dict[str, str]
        { node_name : amino_acid_sequence }
    """
    seqs = defaultdict(str)
    with open(state_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            node, aa = parts[0], parts[2]
            seqs[node] += aa
    return seqs


def parse_indel_file(indel_file):
    """
    Parse the RAxML indel output (already mapped to IQ-TREE node names),
    typically something like:
        Node123  1010100111...

    Returns
    -------
    dict[str, list[str]]
        { node_name : list_of_bits("0"/"1") }
    """
    indels = {}
    with open(indel_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, bits = parts[0], parts[1].strip()
            indels[name] = list(bits)
    return indels


def merge_state_and_indel(state_dict, indel_dict):
    """
    Merge amino-acid ASR (.state) with indel 0/1 information.

    Rules:
      bit = "0" → always '-' (deletion at that column)
      bit = "1" → keep AA as is (if AA is already '-', keep '-')

    Parameters
    ----------
    state_dict : dict[str, str]
        { node_name : AA_sequence }
    indel_dict : dict[str, list[str]]
        { node_name : list_of_bits }

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
                # bit == "1" → keep the AA; if AA is already '-', keep '-'
                new_seq.append(aa if aa != "-" else "-")
        merged[node] = "".join(new_seq)

    if missing_indel:
        print(
            f"[WARN] {len(missing_indel)} nodes found in .state "
            f"but missing from indel file (skipped)."
        )

    # (Optional) nodes that exist in indel_dict but not in state_dict
    extra_indel_nodes = set(indel_dict.keys()) - set(state_dict.keys())
    if extra_indel_nodes:
        print(
            f"[WARN] {len(extra_indel_nodes)} nodes found in indel file "
            f"but not in .state (ignored)."
        )

    return merged


def write_fasta(seqs, outpath, remove_gaps=False):
    """
    Write sequences to FASTA.

    Parameters
    ----------
    seqs : dict[str, str]
        { name : sequence }
    outpath : str
        Output FASTA path.
    remove_gaps : bool, optional
        If True, remove all '-' characters before writing.
    """
    with open(outpath, "w") as f:
        for name, seq in seqs.items():
            out_seq = seq.replace("-", "") if remove_gaps else seq
            f.write(f">{name}\n{out_seq}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine IQ-TREE .state and RAxML-based indel results into final "
            "indel-aware ASR FASTA (with and without gaps)."
        )
    )
    parser.add_argument("--state", required=True,
                        help="IQ-TREE .state file")
    parser.add_argument("--indel", required=True,
                        help="RAxML indel file (already mapped to IQ-TREE node names)")
    parser.add_argument("--out_withgap", default="ASR_final_withgap.fasta",
                        help="Output FASTA (alignment-like, with gaps)")
    parser.add_argument("--out_nogap", default="ASR_final_nogap.fasta",
                        help="Output FASTA (gap-stripped, e.g. AlphaFold input)")
    args = parser.parse_args()

    print(f"[INFO] Reading state file: {args.state}")
    state_dict = parse_state_file(args.state)
    print(f"[INFO] Loaded {len(state_dict)} sequences from state file")

    print(f"[INFO] Reading indel file: {args.indel}")
    indel_dict = parse_indel_file(args.indel)
    print(f"[INFO] Loaded {len(indel_dict)} indel profiles")

    print("[INFO] Merging state and indel information...")
    merged = merge_state_and_indel(state_dict, indel_dict)
    print(f"[INFO] Merged {len(merged)} sequences")

    write_fasta(merged, args.out_withgap, remove_gaps=False)
    write_fasta(merged, args.out_nogap, remove_gaps=True)

    print("[DONE] Output written:")
    print(f"           with gaps : {args.out_withgap}")
    print(f"           no gaps   : {args.out_nogap}")


if __name__ == "__main__":
    main()

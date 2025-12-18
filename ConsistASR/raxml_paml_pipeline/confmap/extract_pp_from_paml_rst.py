#!/usr/bin/env python3
import argparse
import sys
import re

def parse_node_index_from_rst(rst_file, paml_node):
    """
    Scan `node #` lines in the RST file to collect the order of internal nodes,
    and return the 0-based index of the given paml_node in that order.
    """
    order = []
    with open(rst_file) as f:
        for line in f:
            if line.startswith("node #"):
                m = re.match(r"node #(\d+)", line)
                if m:
                    order.append(int(m.group(1)))
    if paml_node not in order:
        sys.exit(f"[ERROR] Node {paml_node} was not found in the RST file.")
    idx = order.index(paml_node)
    return idx, order

def parse_pp_for_node(rst_file, node_index):
    """
    From a PAML .rst file, extract the list of posterior probabilities
    for a given node (node_index) over all alignment columns.

    The function looks for lines following the
      'site   Freq   Data:'
    header, where each site line begins with
      <site_id> <some_int> ...
    and contains a series of '(prob)' entries for each node.

    For each site, it collects all '(... )' floats, and then picks
    the value at position `node_index`.
    """
    pp_vals = []
    with open(rst_file) as f:
        for line in f:
            # Only consider lines that start with "site_number  something"
            if re.match(r"^\s*\d+\s+\d+", line):
                parts = line.strip().split(':', 1)
                if len(parts) < 2:
                    continue
                pp_part = parts[1]
                # Collect all numbers inside parentheses (one per node)
                m = re.findall(r"\(([0-9\.]+)\)", pp_part)
                if m and node_index < len(m):
                    pp_vals.append(float(m[node_index]))
                else:
                    # If parsing fails, append 0.0 as a placeholder
                    pp_vals.append(0.0)
    return pp_vals

def read_withgap_length(path):
    """
    Return the alignment length (number of columns) of the first sequence
    in the with-gap FASTA file.

    This is used only for a length consistency check; indel masks are NOT
    constructed here.
    """
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if seq:
                    break  # only use the first sequence
                continue
            seq.append(line.strip())
    return len(''.join(seq))

def main():
    ap = argparse.ArgumentParser(
        description="Extract per-site PP for a PAML node from .rst (one value per alignment column)."
    )
    ap.add_argument('--rst', required=True, help='PAML .rst file')
    ap.add_argument('--paml_node', type=int, required=True,
                    help='PAML node ID (integer, e.g. 233)')
    ap.add_argument(
        '--withgap_fasta',
        help='with-gap FASTA for this node (used only to check alignment length)'
    )
    ap.add_argument('--out', required=True,
                    help='Output PP file (one value per line)')
    args = ap.parse_args()

    # 1) Determine the node index in the internal-node order
    node_index, order = parse_node_index_from_rst(args.rst, args.paml_node)
    print(f"[INFO] Node {args.paml_node} index: {node_index} / {len(order)}")

    # 2) Extract PP values for this node over all alignment columns
    pp_vals = parse_pp_for_node(args.rst, node_index)
    print(f"[INFO] Extracted PP length (PAML alignment columns) = {len(pp_vals)}")

    # 3) If withgap FASTA is provided, check alignment length
    #    and truncate to the shorter length if necessary.
    if args.withgap_fasta:
        aln_len = read_withgap_length(args.withgap_fasta)
        print(f"[INFO] with-gap alignment length = {aln_len}")

        if aln_len != len(pp_vals):
            print(f"[WARN] PAML PP length ({len(pp_vals)}) and with-gap alignment length ({aln_len}) differ. "
                  f"Truncating to the shorter length.")
            L = min(len(pp_vals), aln_len)
            pp_vals = pp_vals[:L]
            print(f"[INFO] PP length after truncation = {len(pp_vals)}")

    # 4) Write out one PP value per line
    with open(args.out, "w") as out:
        for v in pp_vals:
            out.write(f"{v:.6f}\n")

    print(f"[OK] Wrote PP (one value per alignment column) to {args.out}")

if __name__ == '__main__':
    main()

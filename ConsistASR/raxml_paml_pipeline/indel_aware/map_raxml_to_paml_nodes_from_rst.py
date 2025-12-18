#!/usr/bin/env python3
# ============================================================
# map_raxml_to_paml_nodes_from_rst.py  (tail-normalization version)
#
# This script matches:
#   - the tree with node labels embedded in a PAML .rst file, and
#   - a RAxML_nodeLabelledRootedTree
# by comparing leaf sets, and then rewrites the node labels in
# RAxML_marginalAncestralStates to the corresponding PAML node IDs.
#
# Here we assume that leaf labels differ as:
#   PAML leaf :  112_OG_WP_010903286
#   RAxML leaf:  OG_WP_010903286
# i.e., PAML prefixes the original sequence IDs with "<number>_".
# We therefore apply a fixed "tail-normalization":
#   112_OG_WP_010903286 -> OG_WP_010903286
# ============================================================

import argparse
import re

# ---------- Minimal Newick parser ----------

class Node:
    def __init__(self, name=None):
        # Leaf or internal node name
        self.name = name
        self.children = []        # list[Node]
        self.leafset = None       # frozenset of leaf names


def tokenize_newick(s):
    tokens = []
    buf = []
    for ch in s.strip():
        if ch in '(),:;':
            if buf:
                tokens.append(''.join(buf))
                buf = []
            tokens.append(ch)
        elif ch.isspace():
            if buf:
                tokens.append(''.join(buf))
                buf = []
        else:
            buf.append(ch)
    if buf:
        tokens.append(''.join(buf))
    return tokens


def parse_subtree(tokens, i=0):
    if tokens[i] == '(':
        # Internal node
        i += 1  # skip '('
        children = []
        while True:
            child, i = parse_subtree(tokens, i)
            children.append(child)
            if tokens[i] == ',':
                i += 1
                continue
            elif tokens[i] == ')':
                i += 1
                break
            else:
                raise ValueError(f"Unexpected token after child: {tokens[i]}")

        # Optional node name
        name = None
        if i < len(tokens) and tokens[i] not in [',', ')', ':', ';']:
            name = tokens[i]
            i += 1

        # Optional branch length (ignored)
        if i < len(tokens) and tokens[i] == ':':
            i += 1
            if i < len(tokens):
                i += 1

        node = Node(name)
        node.children = children
        return node, i
    else:
        # Leaf
        name = tokens[i]
        i += 1
        if i < len(tokens) and tokens[i] == ':':
            i += 1
            if i < len(tokens):
                i += 1
        node = Node(name)
        return node, i


def parse_newick_string(s):
    tokens = tokenize_newick(s)
    root, idx = parse_subtree(tokens, 0)
    return root

# ---------- Utility functions ----------

def collect_leafsets(node):
    """
    Recursively compute and attach a leafset (frozenset of leaf names)
    to each node. Returns the set of leaves under this node.
    """
    if not node.children:
        node.leafset = frozenset([node.name])
        return {node.name}
    leaves = set()
    for ch in node.children:
        leaves |= collect_leafsets(ch)
    node.leafset = frozenset(leaves)
    return leaves


def iter_nodes(node):
    """Preorder traversal over all nodes."""
    yield node
    for ch in node.children:
        yield from iter_nodes(ch)


def get_leaf_names(node, acc=None):
    """Collect all leaf names under the given node."""
    if acc is None:
        acc = set()
    if not node.children:
        acc.add(node.name)
    else:
        for ch in node.children:
            get_leaf_names(ch, acc)
    return acc


def transform_leaf_names_tail(root):
    """
    Apply tail-normalization to leaf names:
        '112_OG_WP_010903286' -> 'OG_WP_010903286'
    i.e., drop the leading "<digits>_" prefix when present.
    """
    for nd in iter_nodes(root):
        if nd.children:
            continue
        name = nd.name
        m = re.match(r'^\d+_(.+)$', name)
        if m:
            nd.name = m.group(1)

# ---------- Extract PAML tree from .rst ----------

def extract_paml_tree_from_rst(rst_path):
    """
    Extract the Newick tree (with node labels) from a PAML .rst file.

    We look for the line containing "tree with node labels" (case-insensitive),
    then concatenate subsequent lines until we have a semicolon ';'.
    """
    lines = open(rst_path).read().splitlines()

    start_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "tree with node labels" in low:
            start_idx = i + 1
            break

    if start_idx is None:
        raise RuntimeError(
            f"'tree with node labels' section not found in {rst_path}"
        )

    buf = []
    for line in lines[start_idx:]:
        line = line.strip()
        if not line:
            continue
        buf.append(line)
        if ";" in line:
            break

    if not buf:
        raise RuntimeError("No Newick string found after 'tree with node labels'.")

    tree_str = " ".join(buf)
    if ";" in tree_str:
        tree_str = tree_str.split(";")[0] + ";"

    return tree_str

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Match a RAxML nodeLabelledRootedTree and the PAML tree "
            "embedded in a .rst file by leaf sets (after tail-normalization), "
            "rewrite RAxML_marginalAncestralStates node labels to PAML "
            "node labels, and optionally output a TSV mapping table."
        )
    )
    ap.add_argument("--rst", required=True,
                    help="PAML .rst file containing 'tree with node labels'")
    ap.add_argument("--raxml", required=True,
                    help="RAxML_nodeLabelledRootedTree.* (internal nodes with numeric labels)")
    ap.add_argument("--ancestral", required=True,
                    help="RAxML_marginalAncestralStates.* file (indel states or similar)")
    ap.add_argument("--out", required=True,
                    help="Output file: ancestralStates with PAML node labels")
    ap.add_argument("--map_table", default=None,
                    help="Optional TSV file for RAxML internal node ↔ PAML node mapping")
    ap.add_argument(
        "--paml_tree_out",
        default=None,
        help="Optional output file for the PAML Newick tree with node labels "
             "extracted directly from the .rst file (single-line Newick)."
    )
    args = ap.parse_args()

    # --- PAML tree ---
    paml_tree_str = extract_paml_tree_from_rst(args.rst)

    if args.paml_tree_out:
        with open(args.paml_tree_out, "w") as tf:
            tf.write(paml_tree_str.strip() + "\n")
        print(f"[INFO] Wrote PAML tree with node labels to: {args.paml_tree_out}")

    paml_root = parse_newick_string(paml_tree_str)

    # Leaf names before tail-normalization (for debugging)
    paml_leaf_raw = get_leaf_names(paml_root)
    print(f"[INFO] PAML leaves (raw) example: {sorted(list(paml_leaf_raw))[:5]}")

    # Apply tail-normalization to PAML leaf names
    transform_leaf_names_tail(paml_root)

    # --- RAxML tree ---
    with open(args.raxml) as f:
        rax_str = f.read().strip()
    rax_root = parse_newick_string(rax_str)

    # Collapse dummy ROOT node if present (common in RAxML output)
    if rax_root.name == "ROOT" and len(rax_root.children) == 1:
        print("[INFO] Collapsing ROOT dummy node in RAxML tree")
        rax_root = rax_root.children[0]

    # Collect leaf names and check overlap (for debugging)
    paml_leaf_tail = get_leaf_names(paml_root)
    rax_leaf = get_leaf_names(rax_root)

    print(f"[INFO] PAML leaves (tail) example: {sorted(list(paml_leaf_tail))[:5]}")
    print(f"[INFO] RAxML leaves example:       {sorted(list(rax_leaf))[:5]}")
    print(f"[INFO] Leaf-name overlap size:     {len(paml_leaf_tail & rax_leaf)}")

    # --- Compute leafsets for all nodes ---
    collect_leafsets(paml_root)
    collect_leafsets(rax_root)

    # PAML: leafset -> internal node label
    paml_clade_to_name = {}
    for nd in iter_nodes(paml_root):
        if nd.children and nd.name:
            paml_clade_to_name[nd.leafset] = nd.name

    print(f"[INFO] PAML internal nodes with labels: {len(paml_clade_to_name)}")

    # RAxML: leafset -> internal node label
    rax_clade_to_name = {}
    for nd in iter_nodes(rax_root):
        if nd.children and nd.name:
            rax_clade_to_name[nd.leafset] = nd.name

    print(f"[INFO] RAxML internal nodes with labels: {len(rax_clade_to_name)}")

    # --- Match by leaf set: RAxML → PAML ---
    mapping = {}   # RAxML node label -> PAML node label
    unmapped = []

    for leafset, rname in rax_clade_to_name.items():
        paml_name = paml_clade_to_name.get(leafset)
        if paml_name is None:
            unmapped.append(rname)
        else:
            mapping[rname] = paml_name

    print(f"[INFO] Mapped internal nodes: {len(mapping)}")
    if unmapped:
        print(
            f"[WARN] {len(unmapped)} RAxML nodes could not be mapped "
            f"to any PAML node (topology mismatch or polytomy, etc.)"
        )

    # --- Optional: write mapping table (TSV) ---
    if args.map_table:
        print(f"[INFO] Writing node mapping table to: {args.map_table}")
        with open(args.map_table, "w") as mt:
            mt.write("paml_node\traxml_node\tn_tips\tstatus\n")

            def sort_key(item):
                leafset, rname = item
                if rname.isdigit():
                    return (0, int(rname))
                return (1, rname)

            for leafset, rname in sorted(rax_clade_to_name.items(), key=sort_key):
                paml_name = paml_clade_to_name.get(leafset)
                status = "mapped" if paml_name is not None else "unmapped"
                paml_name_out = paml_name if paml_name is not None else ""
                mt.write(f"{paml_name_out}\t{rname}\t{len(leafset)}\t{status}\n")

    # --- Rewrite RAxML_marginalAncestralStates node labels ---
    n_total = 0
    n_changed = 0

    with open(args.ancestral) as fin, open(args.out, "w") as fout:
        for line in fin:
            stripped = line.strip()
            if not stripped:
                fout.write(line)
                continue

            parts = stripped.split()
            old_name = parts[0]
            new_name = mapping.get(old_name, old_name)
            if new_name != old_name:
                n_changed += 1
            n_total += 1
            fout.write(new_name + "\t" + "\t".join(parts[1:]) + "\n")

    print(f"[DONE] ancestralStates written to: {args.out}")
    print(f"[INFO] Lines processed: {n_total}, renamed: {n_changed}")


if __name__ == "__main__":
    main()

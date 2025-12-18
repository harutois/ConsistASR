#!/usr/bin/env python3
# ============================================================
# map_raxml_to_iqtree_nodes.py
#
# Compare a RAxML nodeLabelledRootedTree and an IQ-TREE
# (already rerooted, with internal nodes labeled as NodeXX),
# match internal nodes that share the same *leaf set*,
# and rewrite RAxML_ancestralStates node labels to IQ-TREE
# NodeXX labels.
#
# Optionally, also write a TSV table describing the mapping.
#
# Example usage:
#   python map_raxml_to_iqtree_nodes.py \
#       --iqtree ASR_QpfamR7_PSI_TM_Coffee.rerooted.nwk \
#       --raxml RAxML_nodeLabelledRootedTree.ASR_QpfamR7_binary \
#       --ancestral RAxML_ancestralStates.ASR_QpfamR7_binary \
#       --out ASR_QpfamR7_PSI_TM_Coffee_indel_ASR.iqtree_named.txt \
#       --map_table ASR_QpfamR7_PSI_TM_Coffee_node_map.tsv
# ============================================================

import argparse
from collections import defaultdict  # (currently unused but harmless)

# ---------- Minimal Newick parser (custom implementation) ----------

class Node:
    def __init__(self, name=None):
        # Leaf or internal node name (SzR_..., HeR_..., "240", "Node12", etc.)
        self.name = name
        self.children = []   # list[Node]
        self.leafset = None  # will be assigned later (frozenset of leaf names)


def normalize_iqtree_label(name: str) -> str:
    """
    Normalize IQ-TREE internal node labels.

    Examples:
      "Node6/100/100" -> "Node6"
      "Node12"        -> "Node12"
      Others          -> As-is
    """
    if name is None:
        return None
    if name.startswith("Node"):
        return name.split('/')[0]
    return name


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
    """
    Recursive Newick parser.
    Returns: (Node, next_index)
    """
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
            i += 1  # skip ':'
            if i < len(tokens):
                i += 1  # skip length token

        node = Node(name)
        node.children = children
        return node, i
    else:
        # Leaf node
        name = tokens[i]
        i += 1
        # Optional branch length (ignored)
        if i < len(tokens) and tokens[i] == ':':
            i += 1
            if i < len(tokens):
                i += 1
        node = Node(name)
        return node, i


def parse_newick_string(s):
    tokens = tokenize_newick(s)
    root, idx = parse_subtree(tokens, 0)
    # trailing ';' (if any) is ignored
    return root


# ---------- Utility functions ----------

def collect_leafsets(node):
    """
    Attach a leafset (frozenset of leaf names) to each node.
    Returns: the set of leaves under this node.
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
    """Traverse the tree and yield all nodes (preorder)."""
    yield node
    for ch in node.children:
        yield from iter_nodes(ch)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Match RAxML nodeLabelledRootedTree and an IQ-TREE rerooted tree "
            "by leaf sets, rewrite RAxML ancestralStates to IQ-TREE NodeXX "
            "labels, and optionally output a TSV mapping table."
        )
    )
    ap.add_argument("--iqtree", required=True,
                    help="IQ-TREE rerooted tree in Newick format (internal nodes labeled as NodeXX)")
    ap.add_argument("--raxml", required=True,
                    help="RAxML_nodeLabelledRootedTree.* (internal nodes with numeric labels)")
    ap.add_argument("--ancestral", required=True,
                    help="RAxML_ancestralStates.* file")
    ap.add_argument("--out", required=True,
                    help="Output file: ancestralStates with IQ-TREE NodeXX labels")
    ap.add_argument("--map_table", default=None,
                    help="Optional TSV file for RAxML internal node ↔ IQ-TREE NodeXX mapping")
    args = ap.parse_args()

    # --- Parse IQ-TREE side ---
    with open(args.iqtree) as f:
        iq_str = f.read().strip()
    iq_root = parse_newick_string(iq_str)
    collect_leafsets(iq_root)

    # IQ-TREE: leafset -> NodeXX name
    iq_clade_to_name = {}
    for nd in iter_nodes(iq_root):
        if nd.children and nd.name:
            norm_name = normalize_iqtree_label(nd.name)
            if norm_name and norm_name.startswith("Node"):
                iq_clade_to_name[nd.leafset] = norm_name

    print(f"[INFO] IQ-TREE internal nodes with NodeXX labels: {len(iq_clade_to_name)}")

    # --- Parse RAxML side ---
    with open(args.raxml) as f:
        rax_str = f.read().strip()
    rax_root = parse_newick_string(rax_str)

    # Collapse dummy ROOT node if it has only one child
    if rax_root.name == "ROOT" and len(rax_root.children) == 1:
        print("[INFO] Collapsing ROOT dummy node in RAxML tree")
        rax_root = rax_root.children[0]

    collect_leafsets(rax_root)

    # RAxML: leafset -> internal node label (numeric)
    rax_clade_to_name = {}
    for nd in iter_nodes(rax_root):
        if (
            nd.children
            and nd.name
            and not nd.name.startswith("SzR_")
            and not nd.name.startswith("HeR_")
            and not nd.name.startswith("OG_")
        ):
            # Internal node label (e.g., "229", "230", ...)
            rax_clade_to_name[nd.leafset] = nd.name

    print(f"[INFO] RAxML internal nodes with numeric labels: {len(rax_clade_to_name)}")

    # --- Match by leaf set (RAxML → IQ-TREE) ---
    mapping = {}   # RAxML node label -> IQ-TREE NodeXX
    unmapped = []  # RAxML node labels that could not be mapped

    for leafset, rname in rax_clade_to_name.items():
        iq_name = iq_clade_to_name.get(leafset)
        if iq_name is None:
            unmapped.append(rname)
        else:
            mapping[rname] = iq_name

    print(f"[INFO] Mapped internal nodes: {len(mapping)}")
    if unmapped:
        print(
            f"[WARN] {len(unmapped)} RAxML nodes could not be mapped "
            f"to any NodeXX (topology mismatch or polytomy near outgroup?)"
        )
        # If needed, you could print the labels:
        # print("[WARN] Unmapped RAxML node labels:", ", ".join(unmapped))

    # --- Optional: write mapping table as TSV ---
    if args.map_table:
        print(f"[INFO] Writing node mapping table to: {args.map_table}")
        with open(args.map_table, "w") as mt:
            mt.write("iqtree_node\traxml_node\tn_tips\tstatus\n")

            # Sort RAxML internal nodes: numeric labels first (by value), then others (by string)
            def sort_key(item):
                leafset, rname = item
                if rname.isdigit():
                    return (0, int(rname))   # numeric first
                return (1, rname)           # then non-numeric, lexicographically

            for leafset, rname in sorted(rax_clade_to_name.items(), key=sort_key):
                iq_name = iq_clade_to_name.get(leafset)
                status = "mapped" if iq_name is not None else "unmapped"
                iq_name_out = iq_name if iq_name is not None else ""
                mt.write(f"{iq_name_out}\t{rname}\t{len(leafset)}\t{status}\n")

    # --- Rewrite ancestralStates ---
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

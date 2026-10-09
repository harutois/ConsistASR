#!/usr/bin/env python3
# ============================================================
# map_raxml_to_iqtree_nodes.py
#
# Compare a RAxML node-labelled rooted tree and an IQ-TREE ASR
# tree (with internal nodes labeled as NodeXX), match internal
# nodes that share the same set of descendant tip labels, and
# rewrite RAxML marginal ancestral-state node labels to IQ-TREE
# NodeXX labels.
#
# Optionally, also write a TSV table describing the mapping.
#
# Example usage:
#   python map_raxml_to_iqtree_nodes.py \
#       --iqtree ASR_example.treefile \
#       --raxml RAxML_nodeLabelledRootedTree.ASR_example_binary \
#       --ancestral RAxML_marginalAncestralStates.ASR_example_binary \
#       --out ASR_example_indel_ASR.iqtree_named.txt \
#       --map_table ASR_example_node_map.tsv
# ============================================================

import argparse
import sys


# ---------- Minimal Newick parser (custom implementation) ----------

class Node:
    def __init__(self, name=None):
        # Leaf or internal node name
        # (e.g. sequence ID, "240", "Node12")
        self.name = name
        self.children = []
        self.leafset = None


def normalize_iqtree_label(name: str) -> str:
    """
    Normalize IQ-TREE internal node labels.

    Examples:
      "Node6/100/100" -> "Node6"
      "Node12"        -> "Node12"
      Others          -> unchanged
    """
    if name is None:
        return None

    if name.startswith("Node"):
        return name.split("/")[0]

    return name


def tokenize_newick(s):
    """Tokenize a simple Newick string."""
    tokens = []
    buf = []

    for ch in s.strip():
        if ch in "(),:;":
            if buf:
                tokens.append("".join(buf))
                buf = []
            tokens.append(ch)

        elif ch.isspace():
            if buf:
                tokens.append("".join(buf))
                buf = []

        else:
            buf.append(ch)

    if buf:
        tokens.append("".join(buf))

    return tokens


def parse_subtree(tokens, i=0):
    """
    Recursive Newick parser.

    Returns
    -------
    node : Node
        Parsed subtree root.
    next_index : int
        Index of the next unconsumed token.
    """
    if i >= len(tokens):
        raise ValueError("Unexpected end of Newick string.")

    if tokens[i] == "(":
        # Internal node
        i += 1
        children = []

        while True:
            child, i = parse_subtree(tokens, i)
            children.append(child)

            if i >= len(tokens):
                raise ValueError(
                    "Unexpected end of Newick string inside internal node."
                )

            if tokens[i] == ",":
                i += 1
                continue

            if tokens[i] == ")":
                i += 1
                break

            raise ValueError(
                f"Unexpected token after child: {tokens[i]}"
            )

        # Optional internal-node name
        name = None
        if (
            i < len(tokens)
            and tokens[i] not in [",", ")", ":", ";"]
        ):
            name = tokens[i]
            i += 1

        # Optional branch length (ignored)
        if i < len(tokens) and tokens[i] == ":":
            i += 1
            if i >= len(tokens):
                raise ValueError(
                    "Missing branch length after ':'."
                )
            i += 1

        node = Node(name)
        node.children = children
        return node, i

    # Leaf node
    name = tokens[i]

    if name in [",", ")", ":", ";", "("]:
        raise ValueError(
            f"Unexpected token where leaf name was expected: {name}"
        )

    i += 1

    # Optional branch length (ignored)
    if i < len(tokens) and tokens[i] == ":":
        i += 1
        if i >= len(tokens):
            raise ValueError(
                "Missing branch length after ':'."
            )
        i += 1

    return Node(name), i


def parse_newick_string(s):
    """Parse a Newick string and return its root Node."""
    tokens = tokenize_newick(s)

    if not tokens:
        raise ValueError("Empty Newick input.")

    root, idx = parse_subtree(tokens, 0)

    # A trailing semicolon is allowed.
    if idx < len(tokens) and tokens[idx] == ";":
        idx += 1

    if idx != len(tokens):
        raise ValueError(
            "Unexpected extra tokens after the Newick tree."
        )

    return root


# ---------- Tree utility functions ----------

def collect_leafsets(node):
    """
    Assign the set of descendant tip labels to every node.

    Returns
    -------
    set
        Leaf names descending from this node.
    """
    if not node.children:
        if node.name is None:
            raise ValueError("Unnamed leaf detected in tree.")

        node.leafset = frozenset([node.name])
        return {node.name}

    leaves = set()

    for child in node.children:
        leaves |= collect_leafsets(child)

    node.leafset = frozenset(leaves)
    return leaves


def iter_nodes(node):
    """Traverse the tree and yield all nodes in preorder."""
    yield node

    for child in node.children:
        yield from iter_nodes(child)


def check_taxon_sets(iq_root, rax_root):
    """
    Confirm that IQ-TREE and RAxML trees contain exactly the same
    set of tip labels.
    """
    iq_tips = iq_root.leafset
    rax_tips = rax_root.leafset

    if iq_tips == rax_tips:
        print(
            f"[INFO] Tree taxon sets match: {len(iq_tips)} taxa"
        )
        return

    only_iq = sorted(iq_tips - rax_tips)
    only_rax = sorted(rax_tips - iq_tips)

    print(
        "[ERROR] IQ-TREE and RAxML trees do not contain "
        "the same taxon set.",
        file=sys.stderr,
    )

    if only_iq:
        print(
            "[ERROR] Taxa found only in IQ-TREE tree: "
            + ", ".join(only_iq),
            file=sys.stderr,
        )

    if only_rax:
        print(
            "[ERROR] Taxa found only in RAxML tree: "
            + ", ".join(only_rax),
            file=sys.stderr,
        )

    sys.exit(1)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Match internal nodes between a RAxML node-labelled "
            "rooted tree and an IQ-TREE ASR tree by descendant "
            "tip sets, rewrite RAxML marginal ancestral-state "
            "node labels to IQ-TREE NodeXX labels, and optionally "
            "write a TSV mapping table."
        )
    )

    ap.add_argument(
        "--iqtree",
        required=True,
        help=(
            "IQ-TREE ASR tree in Newick format "
            "(internal nodes labeled as NodeXX)"
        ),
    )

    ap.add_argument(
        "--raxml",
        required=True,
        help=(
            "RAxML_nodeLabelledRootedTree.* "
            "(internal nodes with numeric labels)"
        ),
    )

    ap.add_argument(
        "--ancestral",
        required=True,
        help="RAxML_marginalAncestralStates.* file",
    )

    ap.add_argument(
        "--out",
        required=True,
        help=(
            "Output file containing RAxML ancestral states "
            "with mapped IQ-TREE NodeXX labels"
        ),
    )

    ap.add_argument(
        "--map_table",
        default=None,
        help=(
            "Optional TSV file containing the RAxML internal "
            "node to IQ-TREE NodeXX mapping"
        ),
    )

    args = ap.parse_args()

    try:
        # ----------------------------------------------------
        # Parse IQ-TREE tree
        # ----------------------------------------------------
        with open(args.iqtree) as f:
            iq_str = f.read().strip()

        iq_root = parse_newick_string(iq_str)
        collect_leafsets(iq_root)

        # IQ-TREE:
        # descendant-tip set -> NodeXX
        iq_clade_to_name = {}

        for nd in iter_nodes(iq_root):
            if nd.children and nd.name:
                norm_name = normalize_iqtree_label(nd.name)

                if norm_name and norm_name.startswith("Node"):
                    iq_clade_to_name[nd.leafset] = norm_name

        print(
            "[INFO] IQ-TREE internal nodes with NodeXX labels: "
            f"{len(iq_clade_to_name)}"
        )

        # ----------------------------------------------------
        # Parse RAxML tree
        # ----------------------------------------------------
        with open(args.raxml) as f:
            rax_str = f.read().strip()

        rax_root = parse_newick_string(rax_str)

        # Some RAxML outputs may contain a dummy ROOT node.
        if (
            rax_root.name == "ROOT"
            and len(rax_root.children) == 1
        ):
            print(
                "[INFO] Collapsing ROOT dummy node in RAxML tree"
            )
            rax_root = rax_root.children[0]

        collect_leafsets(rax_root)
        
        # Read node labels actually present in the
        # RAxML marginal ancestral-state file
        ancestral_labels = set()

        with open(args.ancestral) as fin:
            for line in fin:
                stripped = line.strip()

                if not stripped:
                    continue

                parts = stripped.split()
                if parts:
                    ancestral_labels.add(parts[0])

        print(
            "[INFO] RAxML ancestral-state labels found: "
            f"{len(ancestral_labels)}"
        )

        # ----------------------------------------------------
        # Confirm identical taxon sets before node mapping
        # ----------------------------------------------------
        check_taxon_sets(iq_root, rax_root)

        # RAxML:
        # descendant-tip set -> ancestral internal-node label
        rax_clade_to_name = {}

        for nd in iter_nodes(rax_root):
            if (
                nd.children
                and nd.name
                and nd.name in ancestral_labels
            ):
                rax_clade_to_name[nd.leafset] = nd.name

        print(
            "[INFO] RAxML internal nodes with ancestral states: "
            f"{len(rax_clade_to_name)}"
        )

        # ----------------------------------------------------
        # Match RAxML nodes to IQ-TREE nodes
        # ----------------------------------------------------
        mapping = {}
        unmapped = []

        for leafset, rname in rax_clade_to_name.items():
            iq_name = iq_clade_to_name.get(leafset)

            if iq_name is None:
                unmapped.append(rname)
            else:
                mapping[rname] = iq_name

        print(
            f"[INFO] Mapped internal nodes: {len(mapping)}"
        )

        if unmapped:
            print(
                f"[WARN] {len(unmapped)} RAxML internal nodes "
                "could not be mapped to an IQ-TREE NodeXX label "
                "(for example because of differences in internal-"
                "node resolution near a polytomy or root).",
                file=sys.stderr,
            )

        # ----------------------------------------------------
        # Write optional node-mapping table
        # ----------------------------------------------------
        if args.map_table:
            print(
                "[INFO] Writing node mapping table to: "
                f"{args.map_table}"
            )

            with open(args.map_table, "w") as mt:
                mt.write(
                    "iqtree_node\traxml_node\tn_tips\tstatus\n"
                )

                def sort_key(item):
                    leafset, rname = item

                    if rname.isdigit():
                        return (0, int(rname))

                    return (1, rname)

                for leafset, rname in sorted(
                    rax_clade_to_name.items(),
                    key=sort_key,
                ):
                    iq_name = iq_clade_to_name.get(leafset)

                    status = (
                        "mapped"
                        if iq_name is not None
                        else "unmapped"
                    )

                    iq_name_out = (
                        iq_name
                        if iq_name is not None
                        else ""
                    )

                    mt.write(
                        f"{iq_name_out}\t"
                        f"{rname}\t"
                        f"{len(leafset)}\t"
                        f"{status}\n"
                    )

        # ----------------------------------------------------
        # Rewrite RAxML marginal ancestral-state node labels
        # ----------------------------------------------------
        n_total = 0
        n_changed = 0

        with open(args.ancestral) as fin, open(
            args.out, "w"
        ) as fout:

            for line in fin:
                stripped = line.strip()

                if not stripped:
                    fout.write(line)
                    continue

                parts = stripped.split()

                if not parts:
                    continue

                old_name = parts[0]
                new_name = mapping.get(old_name, old_name)

                if new_name != old_name:
                    n_changed += 1

                n_total += 1

                if len(parts) > 1:
                    fout.write(
                        new_name
                        + "\t"
                        + "\t".join(parts[1:])
                        + "\n"
                    )
                else:
                    fout.write(new_name + "\n")

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

    print(
        f"[DONE] Marginal ancestral states written to: {args.out}"
    )
    print(
        f"[INFO] Lines processed: {n_total}, "
        f"renamed: {n_changed}"
    )


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
# ============================================================
# map_raxml_to_paml_nodes_from_rst.py
#
# Match internal nodes between:
#   - the PAML node-labelled tree embedded in an rst file, and
#   - a RAxML_nodeLabelledRootedTree
#
# Internal nodes are matched by their sets of descendant tip
# labels rather than by node numbers.
#
# PAML may prefix tip labels with "<number>_", for example:
#
#   PAML leaf : 112_OG_WP_010903286
#   RAxML leaf: OG_WP_010903286
#
# Numeric prefixes are removed only when doing so produces an
# exact match to a RAxML tip label.
#
# The script then rewrites node labels in
# RAxML_marginalAncestralStates.* to the corresponding PAML
# node IDs.
# ============================================================

import argparse
import re
import sys


# ---------- Minimal Newick parser ----------

class Node:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.leafset = None


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

    # Leaf
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
    """Parse a Newick string and return its root node."""
    tokens = tokenize_newick(s)

    if not tokens:
        raise ValueError("Empty Newick input.")

    root, idx = parse_subtree(tokens, 0)

    if idx < len(tokens) and tokens[idx] == ";":
        idx += 1

    if idx != len(tokens):
        raise ValueError(
            "Unexpected extra tokens after the Newick tree."
        )

    return root


# ---------- Tree utility functions ----------

def iter_nodes(node):
    """Traverse the tree in preorder."""
    yield node

    for child in node.children:
        yield from iter_nodes(child)


def get_leaf_nodes(node):
    """Return all leaf nodes without discarding duplicate names."""
    return [
        nd
        for nd in iter_nodes(node)
        if not nd.children
    ]


def get_leaf_names(node):
    """Return the set of leaf names in a tree."""
    return {
        nd.name
        for nd in get_leaf_nodes(node)
    }


def collect_leafsets(node):
    """
    Assign the set of descendant tip labels to every node.

    Returns
    -------
    set
        Tip labels descending from this node.
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


# ---------- PAML leaf-name normalization ----------

def normalize_paml_leaf_names(paml_root, raxml_tip_names):
    """
    Normalize PAML leaf names against the actual RAxML tip names.

    Rules
    -----
    1. If the raw PAML name already exactly matches a RAxML tip,
       leave it unchanged.

    2. Otherwise, if the PAML name has a leading "<digits>_"
       prefix and removing that prefix produces an exact RAxML
       tip name, use the prefix-stripped name.

    3. Otherwise, leave the name unchanged. A later taxon-set
       check will report the mismatch.

    Raises
    ------
    ValueError
        If normalization creates duplicate PAML tip names.
    """
    changed = []

    leaf_nodes = get_leaf_nodes(paml_root)

    for nd in leaf_nodes:
        raw_name = nd.name

        if raw_name in raxml_tip_names:
            continue

        match = re.match(r"^\d+_(.+)$", raw_name)

        if match:
            tail = match.group(1)

            if tail in raxml_tip_names:
                nd.name = tail
                changed.append((raw_name, tail))

    normalized_names = [nd.name for nd in leaf_nodes]

    if len(normalized_names) != len(set(normalized_names)):
        duplicates = sorted({
            name
            for name in normalized_names
            if normalized_names.count(name) > 1
        })

        raise ValueError(
            "PAML leaf-name normalization produced duplicate "
            "tip labels: "
            + ", ".join(duplicates)
        )

    print(
        "[INFO] PAML leaf labels normalized by numeric-prefix "
        f"removal: {len(changed)}"
    )

    if changed:
        examples = ", ".join(
            f"{old} -> {new}"
            for old, new in changed[:5]
        )

        print(
            f"[INFO] Normalization example(s): {examples}"
        )


# ---------- Extract PAML tree from rst ----------

def extract_paml_tree_from_rst(rst_path):
    """
    Extract the Newick tree with node labels from a PAML rst file.

    The function searches for the line containing
    'tree with node labels' (case-insensitive), then concatenates
    subsequent non-empty lines until a semicolon is encountered.
    """
    with open(rst_path) as f:
        lines = f.read().splitlines()

    start_idx = None

    for i, line in enumerate(lines):
        if "tree with node labels" in line.lower():
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError(
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
        raise ValueError(
            "No Newick string found after "
            "'tree with node labels'."
        )

    tree_str = " ".join(buf)

    if ";" not in tree_str:
        raise ValueError(
            "PAML node-labelled tree does not contain a "
            "terminating semicolon."
        )

    tree_str = tree_str.split(";", 1)[0] + ";"

    return tree_str


# ---------- Taxon-set check ----------

def check_taxon_sets(paml_root, raxml_root):
    """
    Confirm that normalized PAML and RAxML trees contain exactly
    the same set of tip labels.
    """
    paml_tips = get_leaf_names(paml_root)
    raxml_tips = get_leaf_names(raxml_root)

    if paml_tips == raxml_tips:
        print(
            f"[INFO] Tree taxon sets match: "
            f"{len(paml_tips)} taxa"
        )
        return

    only_paml = sorted(paml_tips - raxml_tips)
    only_raxml = sorted(raxml_tips - paml_tips)

    print(
        "[ERROR] PAML and RAxML trees do not contain "
        "the same taxon set after leaf-name normalization.",
        file=sys.stderr,
    )

    if only_paml:
        print(
            "[ERROR] Taxa found only in PAML tree: "
            + ", ".join(only_paml),
            file=sys.stderr,
        )

    if only_raxml:
        print(
            "[ERROR] Taxa found only in RAxML tree: "
            + ", ".join(only_raxml),
            file=sys.stderr,
        )

    raise ValueError(
        "PAML/RAxML taxon-set mismatch."
    )


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Match internal nodes between the PAML node-labelled "
            "tree embedded in an rst file and a RAxML "
            "node-labelled rooted tree by descendant tip sets, "
            "then rewrite RAxML marginal ancestral-state labels "
            "to PAML node IDs."
        )
    )

    ap.add_argument(
        "--rst",
        required=True,
        help=(
            "PAML rst file containing the "
            "'tree with node labels' section"
        ),
    )

    ap.add_argument(
        "--raxml",
        required=True,
        help="RAxML_nodeLabelledRootedTree.* file",
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
            "Output ancestral-state file with PAML node IDs"
        ),
    )

    ap.add_argument(
        "--map_table",
        default=None,
        help=(
            "Optional TSV file containing the "
            "RAxML-to-PAML node mapping"
        ),
    )

    ap.add_argument(
        "--paml_tree_out",
        default=None,
        help=(
            "Optional output file for the original PAML "
            "node-labelled Newick tree extracted from rst"
        ),
    )

    args = ap.parse_args()

    try:
        # ----------------------------------------------------
        # Extract and parse PAML tree
        # ----------------------------------------------------
        paml_tree_str = extract_paml_tree_from_rst(
            args.rst
        )

        if args.paml_tree_out:
            with open(args.paml_tree_out, "w") as tf:
                tf.write(
                    paml_tree_str.strip() + "\n"
                )

            print(
                "[INFO] Wrote PAML tree with node labels to: "
                f"{args.paml_tree_out}"
            )

        paml_root = parse_newick_string(
            paml_tree_str
        )

        paml_leaf_raw = get_leaf_names(paml_root)

        print(
            "[INFO] PAML leaves (raw) example: "
            f"{sorted(paml_leaf_raw)[:5]}"
        )

        # ----------------------------------------------------
        # Parse RAxML tree
        # ----------------------------------------------------
        with open(args.raxml) as f:
            raxml_tree_str = f.read().strip()

        raxml_root = parse_newick_string(
            raxml_tree_str
        )

        # Some RAxML outputs may contain a dummy ROOT node.
        if (
            raxml_root.name == "ROOT"
            and len(raxml_root.children) == 1
        ):
            print(
                "[INFO] Collapsing ROOT dummy node "
                "in RAxML tree"
            )

            raxml_root = raxml_root.children[0]

        raxml_leaf_names = get_leaf_names(
            raxml_root
        )

        print(
            "[INFO] RAxML leaves example:       "
            f"{sorted(raxml_leaf_names)[:5]}"
        )

        # ----------------------------------------------------
        # Normalize PAML tip labels against RAxML tip labels
        # ----------------------------------------------------
        normalize_paml_leaf_names(
            paml_root,
            raxml_leaf_names,
        )

        paml_leaf_normalized = get_leaf_names(
            paml_root
        )

        print(
            "[INFO] PAML leaves (normalized) example: "
            f"{sorted(paml_leaf_normalized)[:5]}"
        )

        # ----------------------------------------------------
        # Confirm identical taxon sets
        # ----------------------------------------------------
        check_taxon_sets(
            paml_root,
            raxml_root,
        )

        # ----------------------------------------------------
        # Compute descendant-tip sets
        # ----------------------------------------------------
        collect_leafsets(paml_root)
        collect_leafsets(raxml_root)

        # PAML:
        # descendant-tip set -> internal-node label
        paml_clade_to_name = {}

        for nd in iter_nodes(paml_root):
            if nd.children and nd.name:
                paml_clade_to_name[
                    nd.leafset
                ] = nd.name

        print(
            "[INFO] PAML internal nodes with labels: "
            f"{len(paml_clade_to_name)}"
        )

        # ----------------------------------------------------
        # Read RAxML ancestral-state labels
        # ----------------------------------------------------
        ancestral_labels = set()

        with open(args.ancestral) as fin:
            for line in fin:
                stripped = line.strip()

                if not stripped:
                    continue

                parts = stripped.split()

                if parts:
                    ancestral_labels.add(
                        parts[0]
                    )

        if not ancestral_labels:
            raise ValueError(
                "No node labels were found in the "
                "RAxML marginal ancestral-state file."
            )

        print(
            "[INFO] RAxML ancestral-state labels found: "
            f"{len(ancestral_labels)}"
        )

        # RAxML:
        # descendant-tip set -> ancestral internal-node label
        raxml_clade_to_name = {}

        for nd in iter_nodes(raxml_root):
            if (
                nd.children
                and nd.name
                and nd.name in ancestral_labels
            ):
                raxml_clade_to_name[
                    nd.leafset
                ] = nd.name

        print(
            "[INFO] RAxML internal nodes with ancestral states: "
            f"{len(raxml_clade_to_name)}"
        )

        # ----------------------------------------------------
        # Match RAxML nodes to PAML nodes
        # ----------------------------------------------------
        mapping = {}
        unmapped_raxml = []

        for (
            leafset,
            raxml_name,
        ) in raxml_clade_to_name.items():

            paml_name = paml_clade_to_name.get(
                leafset
            )

            if paml_name is None:
                unmapped_raxml.append(
                    raxml_name
                )
            else:
                mapping[
                    raxml_name
                ] = paml_name

        print(
            f"[INFO] Mapped internal nodes: "
            f"{len(mapping)}"
        )

        if unmapped_raxml:
            print(
                f"[WARN] {len(unmapped_raxml)} RAxML "
                "internal node(s) could not be mapped to "
                "a PAML node (for example because of "
                "differences in internal-node resolution "
                "near a polytomy or root).",
                file=sys.stderr,
            )

        # Also report PAML internal nodes that received no
        # corresponding RAxML binary ancestral state.
        mapped_paml_names = set(mapping.values())

        unmapped_paml = sorted(
            set(paml_clade_to_name.values())
            - mapped_paml_names
        )

        if unmapped_paml:
            print(
                f"[WARN] {len(unmapped_paml)} PAML "
                "internal node(s) have no matching RAxML "
                "binary ancestral-state node: "
                + ", ".join(unmapped_paml),
                file=sys.stderr,
            )

        # ----------------------------------------------------
        # Write optional mapping table
        # ----------------------------------------------------
        if args.map_table:
            print(
                "[INFO] Writing node mapping table to: "
                f"{args.map_table}"
            )

            with open(args.map_table, "w") as mt:
                mt.write(
                    "paml_node\t"
                    "raxml_node\t"
                    "n_tips\t"
                    "status\n"
                )

                def sort_key(item):
                    leafset, raxml_name = item

                    if raxml_name.isdigit():
                        return (
                            0,
                            int(raxml_name),
                        )

                    return (
                        1,
                        raxml_name,
                    )

                for (
                    leafset,
                    raxml_name,
                ) in sorted(
                    raxml_clade_to_name.items(),
                    key=sort_key,
                ):

                    paml_name = (
                        paml_clade_to_name.get(
                            leafset
                        )
                    )

                    status = (
                        "mapped"
                        if paml_name is not None
                        else "unmapped"
                    )

                    paml_name_out = (
                        paml_name
                        if paml_name is not None
                        else ""
                    )

                    mt.write(
                        f"{paml_name_out}\t"
                        f"{raxml_name}\t"
                        f"{len(leafset)}\t"
                        f"{status}\n"
                    )

        # ----------------------------------------------------
        # Rewrite RAxML marginal ancestral-state labels
        # ----------------------------------------------------
        n_total = 0
        n_changed = 0

        with open(args.ancestral) as fin, open(
            args.out,
            "w",
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

                if old_name in mapping:
                    new_name = mapping[old_name]
                else:
                    new_name = f"RAxML_UNMAPPED_{old_name}"

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
                    fout.write(
                        new_name + "\n"
                    )

    except (OSError, ValueError) as exc:
        print(
            f"[ERROR] {exc}",
            file=sys.stderr,
        )
        sys.exit(1)

    print(
        "[DONE] Marginal ancestral states written to: "
        f"{args.out}"
    )

    print(
        f"[INFO] Lines processed: {n_total}, "
        f"renamed: {n_changed}"
    )


if __name__ == "__main__":
    main()
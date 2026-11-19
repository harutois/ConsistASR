#!/usr/bin/env python3
"""
map_raxml_to_paml_nodes_from_rst.py

Map RAxML-NG binary-ASR internal nodes to PAML ancestral node IDs.

The PAML tree is rooted and contains n-1 ancestral nodes, including a
degree-2 root. RAxML-NG marginal ASR reconstructs the n-2 genuine
internal vertices of the corresponding unrooted bifurcating tree.
Therefore:

- the PAML root is explicitly excluded from binary indel integration;
- every PAML non-root internal vertex must map one-to-one to exactly one
  RAxML-NG ancestral node;
- mapping is root-independent and uses the unordered three-way tip
  partitions surrounding each unrooted internal vertex.

PAML leaf labels such as "112_OG_WP_010903286" are normalized only when
removing the leading numeric prefix yields an exact RAxML-NG tip label.

Inputs
------
--rst              PAML rst containing "tree with node labels"
--raxml-tree       RAxML-NG *.raxml.ancestralTree
--binary           RAxML-NG two-column binary profile file
--posterior        RAxML-NG binary posterior TSV

Outputs
-------
--out-binary       Binary profiles renamed to PAML node IDs
--out-posterior    Binary posterior TSV renamed to PAML node IDs
--map-table        Mapping table including the explicitly excluded PAML root
--paml-tree-out    Optional extracted PAML node-labelled tree
"""

import argparse
import csv
import re
import sys
from pathlib import Path


class Node:
    def __init__(self, name=None):
        self.name = name
        self.children = []


def tokenize_newick(text):
    tokens, buf = [], []
    for ch in text.strip():
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
    if i >= len(tokens):
        raise ValueError("Unexpected end of Newick string.")

    if tokens[i] == "(":
        i += 1
        children = []

        while True:
            child, i = parse_subtree(tokens, i)
            children.append(child)

            if i >= len(tokens):
                raise ValueError("Unclosed Newick subtree.")

            if tokens[i] == ",":
                i += 1
                continue

            if tokens[i] == ")":
                i += 1
                break

            raise ValueError(f"Unexpected Newick token: {tokens[i]}")

        name = None
        if i < len(tokens) and tokens[i] not in [",", ")", ":", ";"]:
            name = tokens[i]
            i += 1

        if i < len(tokens) and tokens[i] == ":":
            i += 1
            if i >= len(tokens):
                raise ValueError("Missing branch length after ':'.")
            i += 1

        node = Node(name)
        node.children = children
        return node, i

    name = tokens[i]
    if name in [",", ")", ":", ";", "("]:
        raise ValueError(f"Unexpected leaf token: {name}")

    i += 1

    if i < len(tokens) and tokens[i] == ":":
        i += 1
        if i >= len(tokens):
            raise ValueError("Missing branch length after ':'.")
        i += 1

    return Node(name), i


def parse_newick_string(text):
    tokens = tokenize_newick(text)
    if not tokens:
        raise ValueError("Empty Newick input.")

    root, idx = parse_subtree(tokens)

    if idx < len(tokens) and tokens[idx] == ";":
        idx += 1

    if idx != len(tokens):
        raise ValueError("Unexpected extra tokens after Newick tree.")

    return root


def read_newick(path):
    return parse_newick_string(Path(path).read_text().strip())


def iter_nodes(node):
    yield node
    for child in node.children:
        yield from iter_nodes(child)


def get_leaf_nodes(node):
    return [nd for nd in iter_nodes(node) if not nd.children]


def get_leaf_names(node):
    names = [nd.name for nd in get_leaf_nodes(node)]
    if any(name is None for name in names):
        raise ValueError("Unnamed leaf detected.")
    if len(names) != len(set(names)):
        raise ValueError("Duplicate tip labels detected.")
    return set(names)


def extract_paml_tree_from_rst(rst_path):
    with open(rst_path, encoding="utf-8") as f:
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
        stripped = line.strip()
        if not stripped:
            continue
        buf.append(stripped)
        if ";" in stripped:
            break

    if not buf:
        raise ValueError(
            "No Newick string found after 'tree with node labels'."
        )

    tree_str = " ".join(buf)
    if ";" not in tree_str:
        raise ValueError("PAML node-labelled tree has no terminating ';'.")

    return tree_str.split(";", 1)[0] + ";"


def normalize_paml_leaf_names(paml_root, raxml_tip_names):
    changed = []
    leaves = get_leaf_nodes(paml_root)

    for nd in leaves:
        raw = nd.name

        if raw in raxml_tip_names:
            continue

        m = re.match(r"^\d+_(.+)$", raw or "")
        if m and m.group(1) in raxml_tip_names:
            nd.name = m.group(1)
            changed.append((raw, nd.name))

    normalized = [nd.name for nd in leaves]
    if len(normalized) != len(set(normalized)):
        raise ValueError(
            "PAML leaf-name normalization produced duplicate tip labels."
        )

    print(
        f"[INFO] PAML leaf labels normalized by numeric-prefix removal: "
        f"{len(changed)}"
    )


def make_unrooted_graph(root, label, require_paml_root=False):
    """
    Convert rooted Newick representation to an undirected tree.

    For PAML, the labelled degree-2 root is explicitly suppressed and returned.
    For RAxML-NG, an unlabeled or formatting degree-2 root is suppressed if present.
    """
    nodes = list(iter_nodes(root))
    neighbors = {n: set() for n in nodes}

    for parent in nodes:
        for child in parent.children:
            neighbors[parent].add(child)
            neighbors[child].add(parent)

    # Strip unary formatting wrappers.
    while root in neighbors and len(neighbors[root]) == 1 and root.children:
        nxt = next(iter(neighbors[root]))
        neighbors[nxt].remove(root)
        del neighbors[root]
        root = nxt

    excluded_root_name = None

    if require_paml_root:
        if root not in neighbors or len(neighbors[root]) != 2:
            raise ValueError(
                "PAML node-labelled tree is expected to have a degree-2 root."
            )
        if root.name is None or not re.fullmatch(r"\d+", root.name):
            raise ValueError(
                f"PAML root must have a numeric node ID; found {root.name!r}."
            )
        excluded_root_name = root.name

        left, right = tuple(neighbors[root])
        neighbors[left].remove(root)
        neighbors[right].remove(root)
        neighbors[left].add(right)
        neighbors[right].add(left)
        del neighbors[root]

        print(
            f"[INFO] PAML root node excluded from binary integration: "
            f"{excluded_root_name}"
        )

    elif root in neighbors and len(neighbors[root]) == 2:
        # RAxML-NG may serialize an unrooted tree with a degree-2 root.
        left, right = tuple(neighbors[root])
        neighbors[left].remove(root)
        neighbors[right].remove(root)
        neighbors[left].add(right)
        neighbors[right].add(left)
        del neighbors[root]
        print(
            f"[INFO] {label}: suppressed degree-2 Newick root "
            f"({root.name or 'unnamed'})"
        )

    leaves = {}
    for node, adj in neighbors.items():
        if len(adj) == 1:
            if not node.name:
                raise ValueError(f"{label}: unnamed leaf.")
            if node.name in leaves:
                raise ValueError(f"{label}: duplicate tip ID {node.name}.")
            leaves[node.name] = node
        elif len(adj) != 3:
            raise ValueError(
                f"{label}: unsupported unrooted degree {len(adj)} at "
                f"{node.name or '<unlabeled>'}; expected a fully resolved tree."
            )

    return neighbors, leaves, excluded_root_name


def component_tips(neighbors, start, blocked):
    seen = {blocked}
    stack = [start]
    tips = set()

    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)

        if len(neighbors[node]) == 1:
            tips.add(node.name)

        stack.extend(neighbors[node] - seen)

    return frozenset(tips)


def vertex_signatures(neighbors):
    result = {}

    for node, adj in neighbors.items():
        if len(adj) == 1:
            continue

        sig = frozenset(
            component_tips(neighbors, nbr, node)
            for nbr in adj
        )

        if len(sig) != 3 or any(not part for part in sig):
            raise ValueError("Invalid three-way internal-vertex partition.")

        if sig in result:
            raise ValueError("Duplicate three-way internal-vertex signature.")

        result[sig] = node

    return result


def parse_binary_profiles(path):
    profiles = {}

    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip():
                continue

            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Malformed binary-profile line {lineno}: {line.strip()}"
                )

            node, bits = parts
            if node in profiles:
                raise ValueError(f"Duplicate binary node label: {node}")

            if not re.fullmatch(r"[01]+", bits):
                raise ValueError(
                    f"Invalid binary profile for {node}: only 0/1 allowed."
                )

            profiles[node] = bits

    if not profiles:
        raise ValueError("No binary profiles found.")

    lengths = {len(x) for x in profiles.values()}
    if len(lengths) != 1:
        raise ValueError("Binary ancestral profiles have inconsistent lengths.")

    return profiles


def parse_binary_posterior(path):
    required = [
        "node",
        "site",
        "raxmlng_state",
        "p_absent",
        "p_present",
        "map_state",
        "merge_state",
        "status",
    ]

    rows = []
    nodes = set()
    seen = set()

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("Binary posterior TSV has no header.")

        missing = [x for x in required if x not in reader.fieldnames]
        if missing:
            raise ValueError(
                "Binary posterior TSV missing required column(s): "
                + ", ".join(missing)
            )

        for row in reader:
            node = row["node"]
            try:
                site = int(row["site"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid binary posterior site for node {node}: "
                    f"{row['site']}"
                ) from exc

            key = (node, site)
            if key in seen:
                raise ValueError(
                    f"Duplicate binary posterior row: {node} site {site}"
                )
            seen.add(key)
            nodes.add(node)
            rows.append(row)

    if not rows:
        raise ValueError("No binary posterior rows found.")

    return required, rows, nodes


def node_sort_key(name):
    if name.isdigit():
        return (0, int(name))
    m = re.fullmatch(r"Node(\d+)", name)
    if m:
        return (1, int(m.group(1)))
    return (2, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rst", required=True)
    ap.add_argument("--raxml-tree", required=True)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--posterior", required=True)
    ap.add_argument("--out-binary", required=True)
    ap.add_argument("--out-posterior", required=True)
    ap.add_argument("--map-table", required=True)
    ap.add_argument("--paml-tree-out", default=None)
    args = ap.parse_args()

    try:
        paml_tree_str = extract_paml_tree_from_rst(args.rst)

        if args.paml_tree_out:
            Path(args.paml_tree_out).write_text(
                paml_tree_str.strip() + "\n",
                encoding="utf-8",
            )
            print(
                f"[INFO] Wrote PAML tree with node labels to: "
                f"{args.paml_tree_out}"
            )

        paml_root = parse_newick_string(paml_tree_str)
        raxml_root = read_newick(args.raxml_tree)

        raxml_tip_names = get_leaf_names(raxml_root)
        normalize_paml_leaf_names(paml_root, raxml_tip_names)
        paml_tip_names = get_leaf_names(paml_root)

        if paml_tip_names != raxml_tip_names:
            raise ValueError(
                "PAML and RAxML-NG tip sets differ after normalization; "
                f"only_paml={sorted(paml_tip_names - raxml_tip_names)[:10]}, "
                f"only_raxmlng={sorted(raxml_tip_names - paml_tip_names)[:10]}"
            )

        print(f"[INFO] Tree taxon sets match: {len(paml_tip_names)} taxa")

        paml_graph, _, paml_root_id = make_unrooted_graph(
            paml_root,
            "PAML",
            require_paml_root=True,
        )
        raxml_graph, _, _ = make_unrooted_graph(
            raxml_root,
            "RAxML-NG",
            require_paml_root=False,
        )

        paml_sigs = vertex_signatures(paml_graph)
        raxml_sigs = vertex_signatures(raxml_graph)

        expected = len(paml_tip_names) - 2

        if len(paml_sigs) != expected:
            raise ValueError(
                f"PAML non-root internal vertex count is {len(paml_sigs)}, "
                f"expected {expected}."
            )

        if len(raxml_sigs) != expected:
            raise ValueError(
                f"RAxML-NG internal vertex count is {len(raxml_sigs)}, "
                f"expected {expected}."
            )

        paml_by_sig = {}
        paml_labels = set()

        for sig, vertex in paml_sigs.items():
            name = vertex.name
            if not name or not re.fullmatch(r"\d+", name):
                raise ValueError(
                    f"Unlabelled/non-numeric PAML non-root internal vertex: "
                    f"{name!r}"
                )
            if name in paml_labels:
                raise ValueError(f"Duplicate PAML node ID: {name}")
            paml_labels.add(name)
            paml_by_sig[sig] = name

        profiles = parse_binary_profiles(args.binary)
        binary_nodes = set(profiles)

        posterior_header, posterior_rows, posterior_nodes = (
            parse_binary_posterior(args.posterior)
        )

        if binary_nodes != posterior_nodes:
            raise ValueError(
                "Binary profile and posterior node sets differ; "
                f"only_profile={sorted(binary_nodes - posterior_nodes)}, "
                f"only_posterior={sorted(posterior_nodes - binary_nodes)}"
            )

        raxml_by_sig = {}
        raxml_labels = set()

        for sig, vertex in raxml_sigs.items():
            name = vertex.name
            if not name:
                raise ValueError("Unlabelled RAxML-NG internal vertex.")
            if name in raxml_labels:
                raise ValueError(f"Duplicate RAxML-NG node label: {name}")
            raxml_labels.add(name)
            raxml_by_sig[sig] = name

        if raxml_labels != binary_nodes:
            raise ValueError(
                "RAxML-NG ancestralTree internal labels do not match binary "
                "ASR labels; "
                f"only_tree={sorted(raxml_labels - binary_nodes)}, "
                f"only_binary={sorted(binary_nodes - raxml_labels)}"
            )

        if set(paml_by_sig) != set(raxml_by_sig):
            raise ValueError(
                "PAML and RAxML-NG unrooted topologies differ; "
                f"unmatched_paml={len(set(paml_by_sig)-set(raxml_by_sig))}, "
                f"unmatched_raxmlng={len(set(raxml_by_sig)-set(paml_by_sig))}"
            )

        mapping = {
            raxml_by_sig[sig]: paml_by_sig[sig]
            for sig in paml_by_sig
        }

        if (
            len(mapping) != expected
            or len(set(mapping.values())) != expected
        ):
            raise ValueError(
                "RAxML-NG to PAML non-root mapping is not complete one-to-one."
            )

        print(
            f"[INFO] Mapped non-root ancestral nodes: "
            f"{len(mapping)}/{expected} (complete)"
        )
        print(
            f"[INFO] PAML root node {paml_root_id} has no RAxML-NG "
            "binary-ASR counterpart and is excluded from indel-aware outputs."
        )

        # Write renamed two-column binary profiles.
        with open(args.out_binary, "w", encoding="utf-8") as f:
            for raxml_name in sorted(mapping, key=node_sort_key):
                paml_name = mapping[raxml_name]
                f.write(f"{paml_name}\t{profiles[raxml_name]}\n")

        # Write renamed posterior table, retaining the original RAxML-NG label.
        out_fields = [
            "paml_node",
            "raxmlng_node",
            "site",
            "raxmlng_state",
            "p_absent",
            "p_present",
            "map_state",
            "merge_state",
            "status",
        ]

        with open(
            args.out_posterior,
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=out_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in posterior_rows:
                raxml_name = row["node"]
                if raxml_name not in mapping:
                    raise ValueError(
                        f"Posterior row has unmapped RAxML-NG node: "
                        f"{raxml_name}"
                    )

                writer.writerow(
                    {
                        "paml_node": mapping[raxml_name],
                        "raxmlng_node": raxml_name,
                        "site": row["site"],
                        "raxmlng_state": row["raxmlng_state"],
                        "p_absent": row["p_absent"],
                        "p_present": row["p_present"],
                        "map_state": row["map_state"],
                        "merge_state": row["merge_state"],
                        "status": row["status"],
                    }
                )

        # Mapping table.
        with open(args.map_table, "w", encoding="utf-8") as f:
            f.write(
                "paml_node\traxmlng_node\tstatus\tunrooted_partitions\n"
            )

            for sig, paml_name in sorted(
                paml_by_sig.items(),
                key=lambda item: int(item[1]),
            ):
                raxml_name = raxml_by_sig[sig]
                sizes = "/".join(
                    str(x) for x in sorted(len(part) for part in sig)
                )
                f.write(
                    f"{paml_name}\t{raxml_name}\tmapped\t{sizes}\n"
                )

            f.write(
                f"{paml_root_id}\t\texcluded_paml_root\t\n"
            )

        print(
            f"[DONE] PAML-named binary profiles written: {args.out_binary}"
        )
        print(
            f"[DONE] PAML-named binary posterior TSV written: "
            f"{args.out_posterior}"
        )
        print(f"[DONE] Node mapping table written: {args.map_table}")

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

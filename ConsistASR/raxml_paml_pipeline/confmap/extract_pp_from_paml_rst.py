#!/usr/bin/env python3

import argparse
import math
import re
import sys


MARGINAL_MARKER = "Prob of best state at each node, listed by site"


def read_rst(path):
    """Read the complete PAML rst file."""
    try:
        with open(path, encoding="utf-8") as fh:
            return fh.readlines()
    except OSError as exc:
        sys.exit(f"[ERROR] Failed to read RST file {path}: {exc}")


def determine_ancestral_node_order(lines):
    """
    Determine the PAML ancestral-node order from a line such as:

        Nodes 19 to 35 are ancestral

    PAML numbers these internal nodes consecutively.  The resulting
    order is also checked against any later 'node #N' sequence lines
    present in the rst file.
    """

    matches = []

    for lineno, line in enumerate(lines, start=1):
        m = re.match(
            r"^\s*Nodes\s+(\d+)\s+to\s+(\d+)\s+are\s+ancestral\s*$",
            line,
            flags=re.IGNORECASE,
        )

        if m:
            matches.append(
                (lineno, int(m.group(1)), int(m.group(2)))
            )

    if not matches:
        sys.exit(
            "[ERROR] Could not find a line of the form "
            "'Nodes N to M are ancestral' in the RST file."
        )

    if len(matches) > 1:
        desc = ", ".join(
            f"line {lineno}: {first}-{last}"
            for lineno, first, last in matches
        )

        sys.exit(
            "[ERROR] Multiple ancestral-node range declarations "
            f"were found in the RST file: {desc}"
        )

    lineno, first_node, last_node = matches[0]

    if first_node > last_node:
        sys.exit(
            f"[ERROR] Invalid ancestral-node range at line {lineno}: "
            f"{first_node} to {last_node}"
        )

    node_order = list(
        range(first_node, last_node + 1)
    )

    # --------------------------------------------------------
    # Optional cross-check using reconstructed-sequence lines:
    #
    #   node #19   ...
    #   node #20   ...
    #   ...
    #
    # If these lines are present, their node set/order should agree
    # with the declared ancestral-node range.
    # --------------------------------------------------------
    reconstructed_nodes = []

    for line in lines:
        m = re.match(
            r"^\s*node\s*#\s*(\d+)\b",
            line,
            flags=re.IGNORECASE,
        )

        if m:
            reconstructed_nodes.append(
                int(m.group(1))
            )

    if reconstructed_nodes:
        if len(reconstructed_nodes) != len(set(reconstructed_nodes)):
            duplicates = sorted(
                node
                for node in set(reconstructed_nodes)
                if reconstructed_nodes.count(node) > 1
            )

            sys.exit(
                "[ERROR] Duplicate reconstructed node label(s) "
                "found in the RST file: "
                + ", ".join(map(str, duplicates))
            )

        if reconstructed_nodes != node_order:
            sys.exit(
                "[ERROR] Reconstructed 'node #N' sequence order "
                "does not match the ancestral-node range.\n"
                f"Declared order:      {node_order}\n"
                f"Reconstructed order: {reconstructed_nodes}"
            )

    return node_order


def locate_marginal_site_block(lines):
    """
    Locate the first site-data row belonging to:

        Prob of best state at each node, listed by site
    """

    marker_indices = [
        i for i, line in enumerate(lines)
        if MARGINAL_MARKER in line
    ]

    if not marker_indices:
        sys.exit(
            "[ERROR] Marginal reconstruction block not found. "
            f"Expected marker: {MARGINAL_MARKER!r}"
        )

    if len(marker_indices) > 1:
        sys.exit(
            f"[ERROR] Multiple marginal reconstruction blocks "
            f"({len(marker_indices)}) found in the RST file."
        )

    start = marker_indices[0] + 1

    for i in range(start, len(lines)):
        if re.match(
            r"^\s*\d+\s+\d+\s+.*?:",
            lines[i],
        ):
            return i

    sys.exit(
        "[ERROR] No site-data rows were found in the "
        "marginal reconstruction block."
    )


def parse_pp_for_node(
    lines,
    first_site_index,
    node_order,
    paml_node,
):
    """
    Extract one posterior probability per alignment column for
    the requested PAML ancestral node.
    """

    if paml_node not in node_order:
        sys.exit(
            f"[ERROR] Node {paml_node} is not an ancestral node "
            f"in this RST file. "
            f"Available nodes: {node_order[0]}-{node_order[-1]}"
        )

    node_index = node_order.index(paml_node)

    pp_by_site = {}

    prob_re = re.compile(
        r"\("
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
        r"(?:[Ee][-+]?\d+)?)"
        r"\)"
    )

    started = False

    for lineno0 in range(first_site_index, len(lines)):
        line = lines[lineno0]
        lineno = lineno0 + 1

        m = re.match(
            r"^\s*(\d+)\s+(\d+)\s+.*?:(.*)$",
            line,
        )

        if not m:
            if started:
                if not line.strip():
                    continue
                break
            continue

        started = True

        site = int(m.group(1))
        pp_part = m.group(3)

        if site <= 0:
            sys.exit(
                f"[ERROR] Invalid site number at RST line "
                f"{lineno}: {site}"
            )

        if site in pp_by_site:
            sys.exit(
                f"[ERROR] Duplicate site {site} in the "
                "marginal reconstruction block."
            )

        probs_text = prob_re.findall(pp_part)

        if len(probs_text) != len(node_order):
            sys.exit(
                f"[ERROR] Unexpected number of node probabilities "
                f"at site {site} (RST line {lineno}): "
                f"found {len(probs_text)}, expected "
                f"{len(node_order)}."
            )

        try:
            probs = [
                float(value)
                for value in probs_text
            ]
        except ValueError:
            sys.exit(
                f"[ERROR] Failed to parse posterior probabilities "
                f"at site {site}."
            )

        for node_id, value in zip(node_order, probs):
            if not math.isfinite(value):
                sys.exit(
                    f"[ERROR] Non-finite posterior probability "
                    f"for node {node_id} at site {site}: {value}"
                )

            if not (0.0 <= value <= 1.00001):
                sys.exit(
                    f"[ERROR] Posterior probability outside "
                    f"the expected 0-1 range for node {node_id} "
                    f"at site {site}: {value}"
                )

        pp_by_site[site] = probs[node_index]

    if not pp_by_site:
        sys.exit(
            f"[ERROR] No posterior probabilities were extracted "
            f"for node {paml_node}."
        )

    sites = sorted(pp_by_site)

    expected = list(range(1, sites[-1] + 1))

    if sites != expected:
        missing = sorted(
            set(expected) - set(sites)
        )

        preview = ", ".join(
            map(str, missing[:10])
        )

        if len(missing) > 10:
            preview += ", ..."

        sys.exit(
            "[ERROR] Non-contiguous site numbering in the "
            "marginal reconstruction block. "
            f"Missing site(s): {preview}"
        )

    pp_vals = [
        pp_by_site[site]
        for site in sites
    ]

    return node_index, pp_vals


def read_withgap_sequence(path, paml_node):
    """
    Read the relevant aligned ancestral sequence.

    If multiple FASTA entries are present, an exact header
    matching the PAML node ID (e.g. >28) is required.
    Multi-line sequences are supported.
    """

    entries = []
    header = None
    seq_parts = []

    def store():
        if header is not None:
            entries.append(
                (header, "".join(seq_parts))
            )

    try:
        fh = open(path, encoding="utf-8")
    except OSError as exc:
        sys.exit(
            f"[ERROR] Failed to read FASTA file {path}: {exc}"
        )

    with fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")

            if not line.strip():
                continue

            if line.startswith(">"):
                store()
                header = line[1:]
                seq_parts = []
            else:
                if header is None:
                    sys.exit(
                        f"[ERROR] Sequence data found before "
                        f"the first FASTA header at line "
                        f"{lineno} in {path}."
                    )

                seq_parts.append(
                    "".join(line.split())
                )

    store()

    if not entries:
        sys.exit(
            f"[ERROR] No FASTA entries found in {path}."
        )

    target = str(paml_node)

    matches = [
        seq
        for hdr, seq in entries
        if hdr == target
    ]

    if len(matches) > 1:
        sys.exit(
            f"[ERROR] Duplicate exact FASTA header "
            f">{target} found in {path}."
        )

    if len(matches) == 1:
        selected_header = target
        sequence = matches[0]

    elif len(entries) == 1:
        selected_header, sequence = entries[0]

    else:
        sys.exit(
            f"[ERROR] FASTA file {path} contains "
            f"{len(entries)} entries but no exact header "
            f">{target}."
        )

    if not sequence:
        sys.exit(
            f"[ERROR] Empty FASTA sequence for "
            f">{selected_header} in {path}."
        )

    sequence = sequence.upper()

    bad = sorted(
        set(
            c for c in sequence
            if not (c.isalpha() or c == "-")
        )
    )

    if bad:
        sys.exit(
            "[ERROR] Unsupported character(s) in with-gap FASTA: "
            + ", ".join(repr(x) for x in bad)
        )

    return selected_header, sequence


def write_pp(path, pp_vals):
    """Write one posterior probability per line."""

    try:
        with open(path, "w", encoding="utf-8") as out:
            for value in pp_vals:
                out.write(f"{value:.6f}\n")
    except OSError as exc:
        sys.exit(
            f"[ERROR] Failed to write PP file {path}: {exc}"
        )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract site-wise posterior probabilities for one "
            "PAML ancestral node from the marginal reconstruction "
            "block of an rst file."
        )
    )

    ap.add_argument(
        "--rst",
        required=True,
        help="PAML rst file",
    )

    ap.add_argument(
        "--paml_node",
        type=int,
        required=True,
        help="PAML internal-node ID (integer, e.g. 28)",
    )

    ap.add_argument(
        "--withgap_fasta",
        default=None,
        help=(
            "Optional with-gap ancestral FASTA used to verify "
            "alignment-column count. If multiple entries are "
            "present, the exact header >PAML_NODE is selected."
        ),
    )

    ap.add_argument(
        "--out",
        required=True,
        help="Output PP file (one value per alignment column)",
    )

    args = ap.parse_args()

    # --------------------------------------------------------
    # 1) Read RST and determine ancestral-node order
    # --------------------------------------------------------
    lines = read_rst(args.rst)

    node_order = determine_ancestral_node_order(
        lines
    )

    print(
        f"[INFO] PAML ancestral nodes: "
        f"{node_order[0]}-{node_order[-1]} "
        f"(n={len(node_order)})"
    )

    # --------------------------------------------------------
    # 2) Locate marginal reconstruction table
    # --------------------------------------------------------
    first_site_index = locate_marginal_site_block(
        lines
    )

    # --------------------------------------------------------
    # 3) Extract PP for requested node
    # --------------------------------------------------------
    node_index, pp_vals = parse_pp_for_node(
        lines,
        first_site_index,
        node_order,
        args.paml_node,
    )

    print(
        f"[INFO] Node {args.paml_node} index in marginal table: "
        f"{node_index} / {len(node_order)}"
    )

    print(
        f"[INFO] Extracted {len(pp_vals)} PP values "
        f"(sites 1-{len(pp_vals)})."
    )

    # --------------------------------------------------------
    # 4) Optional alignment-coordinate validation
    # --------------------------------------------------------
    if args.withgap_fasta:
        selected_header, sequence = read_withgap_sequence(
            args.withgap_fasta,
            args.paml_node,
        )

        aln_len = len(sequence)

        print(
            f"[INFO] Selected FASTA header: "
            f">{selected_header}"
        )

        print(
            f"[INFO] with-gap alignment length: {aln_len}"
        )

        if aln_len != len(pp_vals):
            sys.exit(
                f"[ERROR] PAML PP length ({len(pp_vals)}) != "
                f"with-gap alignment length ({aln_len}). "
                "The RST output and ancestral FASTA must use "
                "the same alignment-coordinate system."
            )

        print(
            "[INFO] PP/alignment-column counts match."
        )

    # --------------------------------------------------------
    # 5) Write output
    # --------------------------------------------------------
    write_pp(
        args.out,
        pp_vals,
    )

    print(
        f"[OK] Wrote PP (one value per alignment column) "
        f"to {args.out}"
    )


if __name__ == "__main__":
    main()
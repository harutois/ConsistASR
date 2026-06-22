#!/usr/bin/env python3
"""
Score an ESR proxy result from an IQ-TREE .state file.

The reconstructed sequence at an internal proxy node is compared with the
true held-out extant sequence in the original MSA.

Metrics:
- identity over columns where the true target sequence is not a gap
- mean posterior probability assigned to the true residue
"""

import argparse
import sys
from pathlib import Path


def read_fasta_one(path: Path, target: str) -> str:
    name = None
    buf = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name == target:
                    return "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                if name is not None:
                    buf.append(line)

        if name == target:
            return "".join(buf)

    raise SystemExit(f"[ERROR] Target '{target}' not found in {path}")


def parse_iqtree_state(state_path: Path, node_label: str):
    """
    Parse IQ-TREE .state rows for one node.

    Expected columns:
        Node  Site  State  p_A  p_R  ...

    Returns
    -------
    dict[int, tuple[str, dict[str, float]]]
        {site_index_1based: (MAP_state, {AA: posterior_probability})}
    """
    header = None
    rows = {}

    with state_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            if parts[0] == "Node" and parts[1] == "Site" and parts[2] == "State":
                header = parts
                continue

            if header is None:
                raise SystemExit(
                    "[ERROR] Could not find IQ-TREE .state header line "
                    "(expected: Node Site State p_A ...)."
                )

            node = parts[0]
            if node != node_label:
                continue

            try:
                site = int(parts[1])
            except ValueError:
                continue

            state = parts[2]
            probs = {}

            for i in range(3, min(len(header), len(parts))):
                col = header[i]
                if not col.startswith("p_"):
                    continue
                aa = col.replace("p_", "", 1)
                try:
                    probs[aa] = float(parts[i])
                except ValueError:
                    probs[aa] = 0.0

            rows[site] = (state, probs)

    if not rows:
        raise SystemExit(
            f"[ERROR] No rows found for node '{node_label}' in {state_path}. "
            "Check that the ESR proxy tree labels the internal node correctly."
        )

    return rows


def score_esr(true_seq: str, rows, gap_chars: set[str]):
    n_true_nongap = 0
    n_scored = 0
    n_match = 0
    pp_true_sum = 0.0
    missing_sites = 0

    for site, aa_true in enumerate(true_seq, start=1):
        aa_true = aa_true.upper()

        if aa_true in gap_chars:
            continue

        n_true_nongap += 1

        if site not in rows:
            missing_sites += 1
            continue

        aa_pred, probs = rows[site]
        aa_pred = aa_pred.upper()

        n_scored += 1
        if aa_pred == aa_true:
            n_match += 1

        pp_true_sum += probs.get(aa_true, 0.0)

    if n_scored == 0:
        raise SystemExit(
            "[ERROR] No comparable sites were scored. "
            "Check alignment length and .state site numbering."
        )

    identity = 100.0 * n_match / n_scored
    mean_pp_true = pp_true_sum / n_scored

    return {
        "n_true_nongap": n_true_nongap,
        "n_scored": n_scored,
        "n_match": n_match,
        "missing_sites": missing_sites,
        "identity": identity,
        "mean_pp_true": mean_pp_true,
    }


def write_summary(result, target: str, node_label: str, out_handle):
    out_handle.write(f"Target: {target}\n")
    out_handle.write(f"ESR node: {node_label}\n")
    out_handle.write(f"Non-gap sites in true target: {result['n_true_nongap']}\n")
    out_handle.write(f"Comparable sites scored: {result['n_scored']}\n")
    out_handle.write(f"Matching sites: {result['n_match']}\n")
    out_handle.write(f"Sites missing from .state: {result['missing_sites']}\n")
    out_handle.write(f"Identity (%): {result['identity']:.2f}\n")
    out_handle.write(
        f"Mean posterior prob assigned to TRUE residue: {result['mean_pp_true']:.4f}\n"
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Score an ESR proxy result by comparing an IQ-TREE reconstructed "
            "internal node with the original held-out extant sequence."
        )
    )
    parser.add_argument(
        "--orig-msa",
        required=True,
        help="Original MSA FASTA containing the true target sequence.",
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Target extant sequence ID in the original MSA.",
    )
    parser.add_argument(
        "--state",
        required=True,
        help="IQ-TREE .state file from the ESR ASR run.",
    )
    parser.add_argument(
        "--node",
        default=None,
        help="ESR internal node label in the .state file. Default: same as --target.",
    )
    parser.add_argument(
        "--gap-chars",
        default="-?",
        help="Characters treated as gaps/missing in the true sequence. Default: '-?'",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output text file. If omitted, write to stdout.",
    )
    args = parser.parse_args()

    orig_msa = Path(args.orig_msa)
    state_path = Path(args.state)
    target = args.target
    node_label = args.node if args.node is not None else target
    gap_chars = set(args.gap_chars)

    true_seq = read_fasta_one(orig_msa, target)
    rows = parse_iqtree_state(state_path, node_label)

    max_site = max(rows)
    if len(true_seq) != max_site:
        print(
            f"[WARN] Target sequence length ({len(true_seq)}) and maximum .state site "
            f"index ({max_site}) differ. Scoring will use overlapping site indices.",
            file=sys.stderr,
        )

    result = score_esr(true_seq, rows, gap_chars)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as out:
            write_summary(result, target, node_label, out)
    else:
        write_summary(result, target, node_label, sys.stdout)


if __name__ == "__main__":
    main()
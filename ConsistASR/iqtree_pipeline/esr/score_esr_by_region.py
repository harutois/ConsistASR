#!/usr/bin/env python3
"""
Score an ESR proxy result by user-defined sequence regions.

The reconstructed sequence at an ESR internal node is compared with the
true held-out extant sequence in the original MSA, as in score_esr.py.
This script additionally summarizes accuracy by regions defined in an
external TSV file.

Region coordinates are 1-based ungapped residue coordinates of the true
target sequence.

Required region TSV columns:
    region  group  start  end

Example:
    region  group  start  end
    TM1     TM     3      24
    EM1     EM     25     29
"""

import argparse
import csv
import io
import sys
from collections import OrderedDict
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

            if parts[0] != node_label:
                continue

            try:
                site = int(parts[1])
            except ValueError:
                continue

            aa_pred = parts[2]
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

            rows[site] = (aa_pred, probs)

    if not rows:
        raise SystemExit(
            f"[ERROR] No rows found for node '{node_label}' in {state_path}."
        )

    return rows


def read_regions(path: Path):
    """
    Read region definitions from TSV.

    Required columns:
        region, group, start, end

    Coordinates are 1-based ungapped residue coordinates.
    """
    with path.open(encoding="utf-8") as f:
        clean_lines = [
            line for line in f
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not clean_lines:
        raise SystemExit(f"[ERROR] Region file is empty: {path}")

    reader = csv.DictReader(io.StringIO("".join(clean_lines)), delimiter="\t")
    required = {"region", "group", "start", "end"}

    if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
        raise SystemExit(
            "[ERROR] Region TSV must contain columns: region, group, start, end"
        )

    regions = []
    used_positions = {}

    for row in reader:
        name = row["region"].strip()
        group = row["group"].strip()

        try:
            start = int(row["start"])
            end = int(row["end"])
        except ValueError:
            raise SystemExit(f"[ERROR] Invalid start/end in region row: {row}")

        if not name or not group:
            raise SystemExit(f"[ERROR] Empty region/group in row: {row}")

        if start <= 0 or end <= 0 or start > end:
            raise SystemExit(f"[ERROR] Invalid coordinate range in row: {row}")

        for pos in range(start, end + 1):
            if pos in used_positions:
                raise SystemExit(
                    f"[ERROR] Overlapping region definitions at residue {pos}: "
                    f"{used_positions[pos]} and {name}"
                )
            used_positions[pos] = name

        regions.append(
            {
                "region": name,
                "group": group,
                "start": start,
                "end": end,
            }
        )

    return regions


def find_region(residue_index: int, regions):
    for region in regions:
        if region["start"] <= residue_index <= region["end"]:
            return region["region"], region["group"]
    return "UNMAPPED", "UNMAPPED"


def init_stats():
    return {"n": 0, "match": 0, "pp_true": 0.0}


def add_stat(stats, key, match, pp_true):
    if key not in stats:
        stats[key] = init_stats()
    stats[key]["n"] += 1
    stats[key]["match"] += 1 if match else 0
    stats[key]["pp_true"] += pp_true


def format_stat(label, stat):
    n = stat["n"]
    if n == 0:
        return f"{label:<24} sites=0"
    identity = 100.0 * stat["match"] / n
    mean_pp = stat["pp_true"] / n
    return f"{label:<24} sites={n:4d}  identity={identity:6.2f}%  meanPP(true)={mean_pp:.4f}"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Score ESR proxy accuracy by user-defined regions using an "
            "IQ-TREE .state file."
        )
    )
    parser.add_argument("--orig-msa", required=True,
                        help="Original MSA FASTA containing the true target sequence.")
    parser.add_argument("--target", required=True,
                        help="Target extant sequence ID in the original MSA.")
    parser.add_argument("--state", required=True,
                        help="IQ-TREE .state file from the ESR ASR run.")
    parser.add_argument("--regions", required=True,
                        help="Region definition TSV: region, group, start, end.")
    parser.add_argument("--node", default=None,
                        help="ESR internal node label in .state. Default: same as --target.")
    parser.add_argument("--gap-chars", default="-?",
                        help="Characters treated as gaps/missing in the true sequence. Default: '-?'")
    parser.add_argument("--out", default=None,
                        help="Optional output text file. If omitted, write to stdout.")
    args = parser.parse_args()

    target = args.target
    node_label = args.node if args.node is not None else target
    gap_chars = set(args.gap_chars)

    true_seq = read_fasta_one(Path(args.orig_msa), target)
    rows = parse_iqtree_state(Path(args.state), node_label)
    regions = read_regions(Path(args.regions))

    max_site = max(rows)
    if len(true_seq) != max_site:
        print(
            f"[WARN] Target sequence length ({len(true_seq)}) and maximum .state site "
            f"index ({max_site}) differ. Scoring will use overlapping site indices.",
            file=sys.stderr,
        )

    global_stat = init_stats()
    group_stats = OrderedDict()
    region_stats = OrderedDict()
    missing_sites = 0
    non_gap_sites = 0

    for region in regions:
        region_stats[region["region"]] = init_stats()
        if region["group"] not in group_stats:
            group_stats[region["group"]] = init_stats()

    res_idx = 0

    for col_i, aa_true in enumerate(true_seq, start=1):
        aa_true = aa_true.upper()

        if aa_true in gap_chars:
            continue

        non_gap_sites += 1
        res_idx += 1

        if col_i not in rows:
            missing_sites += 1
            continue

        aa_pred, probs = rows[col_i]
        aa_pred = aa_pred.upper()
        match = aa_pred == aa_true
        pp_true = probs.get(aa_true, 0.0)

        region_name, group_name = find_region(res_idx, regions)

        add_stat({"GLOBAL": global_stat}, "GLOBAL", match, pp_true)
        add_stat(group_stats, group_name, match, pp_true)
        add_stat(region_stats, region_name, match, pp_true)

    if global_stat["n"] == 0:
        raise SystemExit("[ERROR] No comparable sites were scored.")

    output = []
    output.append(f"Target: {target}")
    output.append(f"ESR node: {node_label}")
    output.append(f"Region file: {args.regions}")
    output.append(f"Non-gap sites in true target: {non_gap_sites}")
    output.append(f"Comparable sites scored: {global_stat['n']}")
    output.append(f"Sites missing from .state: {missing_sites}")
    output.append("")
    output.append(format_stat("GLOBAL", global_stat))

    output.append("")
    output.append("By group:")
    for group, stat in group_stats.items():
        if stat["n"] > 0:
            output.append(format_stat(group, stat))

    output.append("")
    output.append("By region:")
    for region in regions:
        name = region["region"]
        stat = region_stats[name]
        if stat["n"] > 0:
            output.append(format_stat(name, stat))

    if "UNMAPPED" in region_stats and region_stats["UNMAPPED"]["n"] > 0:
        output.append(format_stat("UNMAPPED", region_stats["UNMAPPED"]))

    text = "\n".join(output) + "\n"

    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
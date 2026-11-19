#!/usr/bin/env python3
"""
raxmlng_binary_asr.py

Parse and validate RAxML-NG binary ancestral-state probabilities for
ConsistASR's shared indel-aware workflows.

Inputs
------
1) RAxML-NG *.raxml.ancestralProbs
2) RAxML-NG *.raxml.ancestralStates
3) Optional IQ-TREE amino-acid ASR .state file

The RAxML-NG probability table is treated as the primary source for
binary presence/absence calls. For each node/site:

    p_1 > p_0  -> present (1)
    p_0 > p_1  -> absent  (0)
    p_0 == p_1 -> ambiguous; WARN and conservatively retain as present (1)

Exact ties are reported explicitly in the posterior TSV and do not stop
the whole workflow. Near-ties are not thresholded: the maximum posterior
state is used as-is, while the probabilities remain available for users
to assess uncertainty.

Validation
----------
Always:
- Node-label sets must match exactly between RAxML-NG ancestralProbs and
  ancestralStates.
- Site sets must match exactly between ancestralProbs and ancestralStates
  for every node.
- Site numbering must be continuous from 1 to N.
- p_0 and p_1 must lie in [0, 1].
- p_0 + p_1 must be approximately 1.
- RAxML-NG ancestralStates and ancestralProbs summary State fields must
  agree exactly.
- When the RAxML-NG summary State is 0 or 1, it must agree with the
  maximum-posterior state. A summary State of "-" is allowed.

When --iqtree-state is supplied:
- RAxML-NG and IQ-TREE node-label sets must match exactly.
- Site sets must match exactly for every node.

Outputs
-------
1) A two-column binary-profile file compatible with downstream FASTA-merging
   scripts:
       Node1<TAB>10101...

2) A per-site TSV with:
       node, site, raxmlng_state, p_absent, p_present,
       map_state, merge_state, status

Here, map_state is 0, 1, or "ambiguous"; merge_state is always 0/1 and
uses 1 for exact ties so that unresolved binary evidence alone does not
delete an amino-acid state.
"""

import argparse
import math
import re
import sys
from collections import defaultdict


PROB_SUM_TOL = 1e-4
TIE_TOL = 1e-12


def natural_node_key(name):
    m = re.fullmatch(r"Node(\d+)", name)
    if m:
        return (0, int(m.group(1)))
    return (1, name)


def parse_iqtree_state(path):
    """Return {node: set(site_numbers)} from an IQ-TREE .state file."""
    sites = defaultdict(set)
    seen_pairs = set()

    with open(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 3:
                continue

            if parts[0] == "Node" and parts[1] == "Site" and parts[2] == "State":
                continue

            node = parts[0]
            try:
                site = int(parts[1])
            except ValueError:
                continue

            key = (node, site)
            if key in seen_pairs:
                raise ValueError(
                    f"IQ-TREE .state contains duplicate node/site row: "
                    f"{node} site {site} (line {line_no})"
                )
            seen_pairs.add(key)
            sites[node].add(site)

    if not sites:
        raise ValueError(f"No ancestral node/site records found in IQ-TREE .state: {path}")

    validate_contiguous_sites(sites, "IQ-TREE .state")
    return dict(sites)


def parse_raxmlng_states(path):
    """Return {node: {site: state_char}} from *.raxml.ancestralStates."""
    states = {}
    with open(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue

            parts = stripped.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Unexpected ancestralStates format at line {line_no}: {stripped}"
                )

            node, seq = parts
            if node in states:
                raise ValueError(f"Duplicate node in ancestralStates: {node}")

            bad = sorted(set(seq) - {"0", "1", "-"})
            if bad:
                raise ValueError(
                    f"Unexpected binary state character(s) for {node}: {bad}"
                )
            states[node] = seq

    if not states:
        raise ValueError(f"No ancestral sequences found in: {path}")

    return states


def parse_raxmlng_probs(path):
    """
    Parse RAxML-NG *.raxml.ancestralProbs.

    Return:
        records[node][site] = {
            "summary_state": str,
            "p0": float,
            "p1": float
        }
    """
    records = defaultdict(dict)

    with open(path) as handle:
        header = None
        column_index = None

        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue

            parts = stripped.split()

            if header is None:
                header = parts
                required = ["Node", "Site", "State", "p_0", "p_1"]
                missing = [x for x in required if x not in header]
                if missing:
                    raise ValueError(
                        f"ancestralProbs header is missing required column(s): "
                        + ", ".join(missing)
                    )
                column_index = {name: header.index(name) for name in required}
                continue

            try:
                node = parts[column_index["Node"]]
                site = int(parts[column_index["Site"]])
                summary_state = parts[column_index["State"]]
                p0 = float(parts[column_index["p_0"]])
                p1 = float(parts[column_index["p_1"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"Could not parse ancestralProbs line {line_no}: {stripped}"
                ) from exc

            if summary_state not in {"0", "1", "-"}:
                raise ValueError(
                    f"Unexpected RAxML-NG State '{summary_state}' for "
                    f"{node} site {site}"
                )

            if site in records[node]:
                raise ValueError(
                    f"Duplicate ancestralProbs row: {node} site {site}"
                )

            for label, value in (("p_0", p0), ("p_1", p1)):
                if not math.isfinite(value) or value < 0.0 or value > 1.0:
                    raise ValueError(
                        f"{label} out of range for {node} site {site}: {value}"
                    )

            if abs((p0 + p1) - 1.0) > PROB_SUM_TOL:
                raise ValueError(
                    f"Posterior probabilities do not sum to ~1 for "
                    f"{node} site {site}: p_0={p0}, p_1={p1}"
                )

            records[node][site] = {
                "summary_state": summary_state,
                "p0": p0,
                "p1": p1,
            }

    if not records:
        raise ValueError(f"No ancestral probability records found in: {path}")

    site_sets = {node: set(site_map) for node, site_map in records.items()}
    validate_contiguous_sites(site_sets, "RAxML-NG ancestralProbs")
    return dict(records)


def validate_contiguous_sites(site_sets, label):
    """Require each node to have sites 1..N without gaps."""
    for node, sites in site_sets.items():
        if not sites:
            raise ValueError(f"{label}: no sites found for {node}")
        expected = set(range(1, max(sites) + 1))
        if sites != expected:
            missing = sorted(expected - sites)
            extra = sorted(sites - expected)
            raise ValueError(
                f"{label}: non-contiguous sites for {node}; "
                f"missing={missing[:10]}, extra={extra[:10]}"
            )


def validate_raxmlng_node_and_site_sets(prob_records, raxml_states):
    """Validate RAxML-NG ancestralProbs against ancestralStates."""
    prob_nodes = set(prob_records)
    state_nodes = set(raxml_states)

    if prob_nodes != state_nodes:
        raise ValueError(
            "RAxML-NG ancestralProbs and ancestralStates node-label sets differ; "
            f"only_probs={sorted(prob_nodes - state_nodes, key=natural_node_key)}, "
            f"only_states={sorted(state_nodes - prob_nodes, key=natural_node_key)}"
        )

    for node in sorted(prob_nodes, key=natural_node_key):
        prob_sites = set(prob_records[node])
        state_sites = set(range(1, len(raxml_states[node]) + 1))

        if prob_sites != state_sites:
            raise ValueError(
                f"RAxML-NG ancestralProbs/ancestralStates site mismatch for {node}: "
                f"probs={len(prob_sites)}, states={len(state_sites)}"
            )


def validate_iqtree_node_and_site_sets(prob_records, iqtree_sites):
    """Validate RAxML-NG binary-ASR coordinates against an IQ-TREE .state file."""
    prob_nodes = set(prob_records)
    iq_nodes = set(iqtree_sites)

    if prob_nodes != iq_nodes:
        raise ValueError(
            "RAxML-NG binary ASR and IQ-TREE .state node-label sets differ; "
            f"only_raxmlng={sorted(prob_nodes - iq_nodes, key=natural_node_key)}, "
            f"only_iqtree={sorted(iq_nodes - prob_nodes, key=natural_node_key)}"
        )

    for node in sorted(prob_nodes, key=natural_node_key):
        prob_sites = set(prob_records[node])
        iq_sites = iqtree_sites[node]

        if prob_sites != iq_sites:
            raise ValueError(
                f"RAxML-NG/IQ-TREE site mismatch for {node}: "
                f"raxmlng={len(prob_sites)}, iqtree={len(iq_sites)}"
            )



def derive_calls(prob_records, raxml_states):
    """
    Return rows and binary profiles.

    For exact ties, map_state='ambiguous' but merge_state='1'.
    """
    rows = []
    profiles = {}
    tie_sites = []
    dash_summary_sites = []

    for node in sorted(prob_records, key=natural_node_key):
        bits = []
        seq = raxml_states[node]

        for site in sorted(prob_records[node]):
            rec = prob_records[node][site]
            p0 = rec["p0"]
            p1 = rec["p1"]
            summary_state = rec["summary_state"]
            state_file_char = seq[site - 1]

            if state_file_char != summary_state:
                raise ValueError(
                    f"RAxML-NG ancestralStates/ancestralProbs State mismatch at "
                    f"{node} site {site}: states='{state_file_char}', "
                    f"probs='{summary_state}'"
                )

            if abs(p1 - p0) <= TIE_TOL:
                map_state = "ambiguous"
                merge_state = "1"
                status = "tie_retained_as_present"
                tie_sites.append((node, site, p0, p1))
            elif p1 > p0:
                map_state = "1"
                merge_state = "1"
                status = "present"
            else:
                map_state = "0"
                merge_state = "0"
                status = "absent"

            if summary_state in {"0", "1"} and summary_state != map_state:
                raise ValueError(
                    f"RAxML-NG summary State disagrees with maximum posterior "
                    f"state at {node} site {site}: State={summary_state}, "
                    f"p_0={p0}, p_1={p1}"
                )

            if summary_state == "-":
                dash_summary_sites.append((node, site, p0, p1))

            bits.append(merge_state)
            rows.append(
                (node, site, summary_state, p0, p1, map_state, merge_state, status)
            )

        profiles[node] = "".join(bits)

    return rows, profiles, tie_sites, dash_summary_sites


def write_binary_profiles(profiles, path):
    with open(path, "w") as handle:
        for node in sorted(profiles, key=natural_node_key):
            handle.write(f"{node}\t{profiles[node]}\n")


def write_posterior_tsv(rows, path):
    with open(path, "w") as handle:
        handle.write(
            "node\tsite\traxmlng_state\tp_absent\tp_present\t"
            "map_state\tmerge_state\tstatus\n"
        )
        for node, site, summary_state, p0, p1, map_state, merge_state, status in rows:
            handle.write(
                f"{node}\t{site}\t{summary_state}\t"
                f"{p0:.5f}\t{p1:.5f}\t{map_state}\t{merge_state}\t{status}\n"
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--probs",
        required=True,
        help="RAxML-NG *.raxml.ancestralProbs file",
    )
    ap.add_argument(
        "--states",
        required=True,
        help="RAxML-NG *.raxml.ancestralStates file",
    )
    ap.add_argument(
        "--iqtree-state",
        default=None,
        help=(
            "Optional IQ-TREE amino-acid ASR .state file. When supplied, "
            "node labels and site coordinates are validated against the "
            "RAxML-NG binary ASR output."
        ),
    )
    ap.add_argument(
        "--out-binary",
        required=True,
        help="Output two-column NodeXX/binary-profile file for FASTA merging",
    )
    ap.add_argument(
        "--out-tsv",
        required=True,
        help="Output per-site posterior probability TSV",
    )
    args = ap.parse_args()

    try:
        prob_records = parse_raxmlng_probs(args.probs)
        raxml_states = parse_raxmlng_states(args.states)

        validate_raxmlng_node_and_site_sets(
            prob_records,
            raxml_states,
        )

        iqtree_sites = None
        if args.iqtree_state is not None:
            iqtree_sites = parse_iqtree_state(args.iqtree_state)
            validate_iqtree_node_and_site_sets(
                prob_records,
                iqtree_sites,
            )

        rows, profiles, tie_sites, dash_summary_sites = derive_calls(
            prob_records,
            raxml_states,
        )

        write_binary_profiles(profiles, args.out_binary)
        write_posterior_tsv(rows, args.out_tsv)

        nodes = sorted(profiles, key=natural_node_key)
        n_sites = len(profiles[nodes[0]])
        total_rows = len(rows)

        print(
            f"[INFO] RAxML-NG ancestralProbs / ancestralStates "
            f"node labels match exactly: {len(nodes)} node(s)"
        )
        if iqtree_sites is not None:
            print(
                f"[INFO] IQ-TREE / RAxML-NG node labels and site coordinates "
                f"match exactly: {len(nodes)} node(s)"
            )
        else:
            print("[INFO] IQ-TREE .state validation: not requested")
        print(f"[INFO] Sites per node: {n_sites}")
        print(f"[INFO] Binary posterior rows validated: {total_rows}")

        if dash_summary_sites:
            print(
                f"[INFO] RAxML-NG summary State='-' at "
                f"{len(dash_summary_sites)} site(s); final binary calls were "
                f"derived from p_0/p_1."
            )

        if tie_sites:
            print(
                f"[WARN] {len(tie_sites)} binary site(s) have equal posterior "
                f"probabilities (p_0 = p_1).",
                file=sys.stderr,
            )
            print(
                "[WARN] These sites are marked 'ambiguous' in the TSV and "
                "retained as present (merge_state=1).",
                file=sys.stderr,
            )
            for node, site, p0, p1 in tie_sites[:10]:
                print(
                    f"[WARN]   {node} site {site}: p_0={p0:.5f}, p_1={p1:.5f}",
                    file=sys.stderr,
                )
            if len(tie_sites) > 10:
                print(
                    f"[WARN]   ... and {len(tie_sites) - 10} more tie site(s)",
                    file=sys.stderr,
                )

        print(f"[DONE] Binary profiles written: {args.out_binary}")
        print(f"[DONE] Binary posterior TSV written: {args.out_tsv}")

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

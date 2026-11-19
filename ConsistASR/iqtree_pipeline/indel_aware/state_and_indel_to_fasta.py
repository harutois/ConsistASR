#!/usr/bin/env python3
# ============================================================
# state_and_indel_to_fasta.py
#
# Combine IQ-TREE .state output and RAxML-NG-derived binary
# gap-state results into indel-aware ancestral FASTA sequences
# and posterior-probability tables.
#
# IQ-TREE may report State="-" when no amino-acid state reaches
# its reporting threshold. This "-" is not interpreted as a
# biological gap. If the binary reconstruction infers presence
# at such a site, the amino acid with the highest posterior
# probability in the IQ-TREE .state row is used in the final
# indel-aware FASTA.
#
# Outputs:
#   1) Raw IQ-TREE ancestral FASTA at alignment length
#      (State column preserved exactly, including "-")
#   2) Raw gap-stripped IQ-TREE ancestral FASTA
#   3) Indel-aware alignment-like FASTA with gaps
#   4) Indel-aware gap-stripped FASTA
#   5) Alignment-site posterior TSV (optional)
#   6) Final nogap residue posterior TSV (optional)
# ============================================================

import argparse
import csv
import math
import sys
from collections import defaultdict


AA_POSTERIOR_PREFIX = "p_"
POSTERIOR_TIE_TOL = 1e-12
POSTERIOR_SUM_TOL = 1e-3


# ---------------------------
# IQ-TREE .state parser
# ---------------------------

def parse_state_file(state_file):
    """
    Parse an IQ-TREE .state file.

    Expected columns
    ----------------
        Node  Site  State  p_A  p_R ...

    Returns
    -------
    raw_seqs : dict[str, str]
        Raw IQ-TREE State-column sequences for FASTA output.

    state_records : dict[str, list[dict]]
        Per-site records in alignment order.

    posterior_labels : list[str]
        State labels corresponding to posterior columns, in header order.

    Notes
    -----
    - Site numbers must be contiguous and 1-based for every node.
    - Raw State values are preserved exactly for raw FASTA output.
    - Posterior columns are identified dynamically from header names
      beginning with "p_".
    """
    header = None
    posterior_columns = []
    posterior_labels = []
    states = defaultdict(list)

    with open(state_file) as f:
        for lineno, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split()

            if (
                len(parts) >= 3
                and parts[0] == "Node"
                and parts[1] == "Site"
                and parts[2] == "State"
            ):
                if header is not None:
                    raise ValueError(
                        "Multiple IQ-TREE .state header lines were detected."
                    )

                header = parts
                posterior_columns = [
                    (idx, col[2:], col)
                    for idx, col in enumerate(header)
                    if col.startswith(AA_POSTERIOR_PREFIX) and len(col) > 2
                ]
                posterior_labels = [label for _, label, _ in posterior_columns]

                if not posterior_columns:
                    raise ValueError(
                        "No posterior-probability columns (p_A, p_R, ...) "
                        "were found in the IQ-TREE .state header."
                    )
                continue

            if header is None:
                raise ValueError(
                    "IQ-TREE .state header was not found before data rows."
                )

            if len(parts) < 3:
                raise ValueError(
                    f"Malformed .state row at line {lineno}: '{line.strip()}'"
                )

            node = parts[0]

            try:
                site = int(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid site number in .state file at line {lineno}: "
                    f"'{parts[1]}'"
                ) from exc

            state = parts[2]

            posteriors = {}
            raw_posteriors = {}

            for idx, label, colname in posterior_columns:
                if idx >= len(parts):
                    raise ValueError(
                        f"Missing posterior value for {colname} at "
                        f"{node} site {site} (line {lineno})."
                    )

                raw_value = parts[idx]

                try:
                    value = float(raw_value)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid posterior value for {colname} at "
                        f"{node} site {site}: '{raw_value}'"
                    ) from exc

                if not math.isfinite(value) or value < 0.0 or value > 1.0:
                    raise ValueError(
                        f"Posterior probability out of range for {colname} at "
                        f"{node} site {site}: {value}"
                    )

                posteriors[label] = value
                raw_posteriors[colname] = raw_value

            total_pp = sum(posteriors.values())
            if abs(total_pp - 1.0) > POSTERIOR_SUM_TOL:
                raise ValueError(
                    f"Amino-acid posterior probabilities do not sum to ~1 at "
                    f"{node} site {site}: sum={total_pp:.8f}"
                )

            states[node].append(
                {
                    "site": site,
                    "state": state,
                    "posteriors": posteriors,
                    "raw_posteriors": raw_posteriors,
                }
            )

    if not states:
        raise ValueError(
            "No ancestral states were found in the IQ-TREE .state file."
        )

    raw_seqs = {}
    state_records = {}

    for node, records in states.items():
        records.sort(key=lambda x: x["site"])

        sites = [rec["site"] for rec in records]
        expected_sites = list(range(1, len(records) + 1))

        if sites != expected_sites:
            raise ValueError(
                f"Non-contiguous or duplicated site numbers detected "
                f"for node {node}. Expected sites 1-{len(records)}, "
                f"but found an inconsistent site series."
            )

        raw_seqs[node] = "".join(rec["state"] for rec in records)
        state_records[node] = records

    return raw_seqs, state_records, posterior_labels


# ---------------------------
# Binary gap-state parser
# ---------------------------

def parse_indel_file(indel_file):
    """
    Parse the two-column binary gap-state file.

    Expected format
    ---------------
        Node123  1010100111...

    Returns
    -------
    dict[str, list[str]]
        {node_name: list_of_bits}
    """
    indels = {}

    with open(indel_file) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 2:
                raise ValueError(
                    f"Malformed binary gap-state line at line {lineno}: '{line}'"
                )

            name = parts[0]
            bits = parts[1].strip()

            if name in indels:
                raise ValueError(
                    f"Duplicate node name in binary gap-state file: '{name}'"
                )

            invalid_chars = set(bits) - {"0", "1"}
            if invalid_chars:
                invalid_display = "".join(sorted(invalid_chars))
                raise ValueError(
                    f"Invalid binary state for node {name}: "
                    f"unexpected character(s) '{invalid_display}'. "
                    "Only '0' and '1' are allowed."
                )

            indels[name] = list(bits)

    if not indels:
        raise ValueError(
            "No binary gap-state profiles were found in the indel file."
        )

    return indels


# ---------------------------
# Binary posterior TSV parser
# ---------------------------

def parse_binary_posterior(binary_posterior_file):
    """
    Parse the TSV written by raxmlng_binary_asr.py.

    Required columns
    ----------------
        node
        site
        raxmlng_state
        p_absent
        p_present
        map_state
        merge_state
        status

    Returns
    -------
    dict[str, dict[int, dict]]
        {node: {site: record}}
    """
    required = {
        "node",
        "site",
        "raxmlng_state",
        "p_absent",
        "p_present",
        "map_state",
        "merge_state",
        "status",
    }

    data = defaultdict(dict)

    with open(binary_posterior_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(
                f"Binary posterior TSV has no header: {binary_posterior_file}"
            )

        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Binary posterior TSV is missing required column(s): "
                + ", ".join(sorted(missing))
            )

        for lineno, row in enumerate(reader, start=2):
            node = row["node"]

            try:
                site = int(row["site"])
                p_absent = float(row["p_absent"])
                p_present = float(row["p_present"])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed binary posterior row at line {lineno}."
                ) from exc

            if site in data[node]:
                raise ValueError(
                    f"Duplicate binary posterior row: {node} site {site}"
                )

            if not (0.0 <= p_absent <= 1.0 and 0.0 <= p_present <= 1.0):
                raise ValueError(
                    f"Binary posterior probability out of range at "
                    f"{node} site {site}."
                )

            merge_state = row["merge_state"]
            if merge_state not in {"0", "1"}:
                raise ValueError(
                    f"Invalid merge_state '{merge_state}' at {node} site {site}."
                )

            data[node][site] = {
                "raxmlng_state": row["raxmlng_state"],
                "p_absent": p_absent,
                "p_present": p_present,
                "map_state": row["map_state"],
                "merge_state": merge_state,
                "status": row["status"],
            }

    if not data:
        raise ValueError(
            "No rows were found in the binary posterior TSV."
        )

    return dict(data)


# ---------------------------
# Amino-acid posterior helpers
# ---------------------------

def get_max_posterior_state(record, posterior_labels):
    """
    Return maximum-posterior state, PP, and all tied states.

    Ties are resolved deterministically by posterior-column order.
    """
    posteriors = record["posteriors"]
    max_pp = max(posteriors.values())

    tied = [
        aa
        for aa in posterior_labels
        if abs(posteriors[aa] - max_pp) <= POSTERIOR_TIE_TOL
    ]

    if not tied:
        raise ValueError(
            f"No maximum-posterior state found at site {record['site']}."
        )

    return tied[0], max_pp, tied


def resolve_undetermined_state(node, record, posterior_labels):
    """
    Resolve IQ-TREE State="-" using the maximum posterior state.
    """
    chosen_aa, chosen_pp, tied = get_max_posterior_state(
        record,
        posterior_labels,
    )

    return chosen_aa, chosen_pp, tied


# ---------------------------
# Consistency checks
# ---------------------------

def validate_node_sets(state_records, indel_dict):
    state_nodes = set(state_records)
    indel_nodes = set(indel_dict)

    if state_nodes != indel_nodes:
        raise ValueError(
            "IQ-TREE .state and binary gap-state node sets differ; "
            f"only_state={sorted(state_nodes - indel_nodes)}, "
            f"only_binary={sorted(indel_nodes - state_nodes)}"
        )


def validate_binary_posterior_against_profiles(
    state_records,
    indel_dict,
    binary_posterior,
):
    state_nodes = set(state_records)
    posterior_nodes = set(binary_posterior)

    if state_nodes != posterior_nodes:
        raise ValueError(
            "IQ-TREE .state and binary posterior node sets differ; "
            f"only_state={sorted(state_nodes - posterior_nodes)}, "
            f"only_binary_posterior={sorted(posterior_nodes - state_nodes)}"
        )

    for node, records in state_records.items():
        bits = indel_dict[node]

        if len(bits) != len(records):
            raise ValueError(
                f"Length mismatch for node {node}: "
                f"state={len(records)}, binary={len(bits)}"
            )

        posterior_sites = set(binary_posterior[node])
        expected_sites = set(range(1, len(records) + 1))

        if posterior_sites != expected_sites:
            raise ValueError(
                f"Binary posterior site set mismatch for {node}."
            )

        for record, bit in zip(records, bits):
            site = record["site"]
            posterior_bit = binary_posterior[node][site]["merge_state"]

            if bit != posterior_bit:
                raise ValueError(
                    f"Binary state mismatch at {node} site {site}: "
                    f"profile={bit}, posterior_tsv={posterior_bit}"
                )


# ---------------------------
# Merge AA and binary states
# ---------------------------

def merge_state_and_indel(state_records, indel_dict, posterior_labels):
    """
    Merge amino-acid states with binary gap-state information.

    Rules
    -----
    binary = 0
        Final state is '-'.

    binary = 1 and IQ-TREE State != '-'
        Retain IQ-TREE State.

    binary = 1 and IQ-TREE State == '-'
        Use maximum-posterior amino-acid state.
    """
    validate_node_sets(state_records, indel_dict)

    merged = {}
    resolution_events = []

    for node, records in state_records.items():
        bits = indel_dict[node]

        if len(bits) != len(records):
            raise ValueError(
                f"Length mismatch for node {node}: "
                f"state={len(records)}, binary={len(bits)}"
            )

        new_seq = []

        for record, bit in zip(records, bits):
            state = record["state"]
            site = record["site"]

            if bit == "0":
                new_seq.append("-")

            elif bit == "1":
                if state != "-":
                    if state not in record["posteriors"]:
                        raise ValueError(
                            f"IQ-TREE State '{state}' has no matching posterior "
                            f"column at {node} site {site}."
                        )
                    new_seq.append(state)
                else:
                    chosen_aa, chosen_pp, tied_aas = resolve_undetermined_state(
                        node,
                        record,
                        posterior_labels,
                    )
                    new_seq.append(chosen_aa)

                    resolution_events.append(
                        {
                            "node": node,
                            "site": site,
                            "chosen_aa": chosen_aa,
                            "chosen_pp": chosen_pp,
                            "tied_aas": tied_aas,
                        }
                    )

            else:
                raise ValueError(
                    f"Unexpected binary state '{bit}' for node {node}."
                )

        merged[node] = "".join(new_seq)

    return merged, resolution_events


# ---------------------------
# Posterior TSV writers
# ---------------------------

def write_site_posterior_tsv(
    state_records,
    indel_dict,
    binary_posterior,
    posterior_labels,
    outpath,
):
    """
    Write alignment-site-level posterior information for every node/site.

    Gap-classified sites are retained. This is the complete audit table.
    """
    aa_columns = [f"p_{aa}" for aa in posterior_labels]

    fieldnames = [
        "node",
        "alignment_site",
        "iqtree_state",
        "aa_map",
        "aa_map_pp",
        "aa_map_tie",
        "p_absent",
        "p_present",
        "binary_state",
        "binary_status",
        "final_state",
        "final_residue",
        "final_aa_pp",
        "aa_source",
    ] + aa_columns

    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for node, records in state_records.items():
            residue_number = 0
            bits = indel_dict[node]

            for record, bit in zip(records, bits):
                site = record["site"]
                state = record["state"]
                b = binary_posterior[node][site]

                aa_map, aa_map_pp, tied = get_max_posterior_state(
                    record,
                    posterior_labels,
                )

                if bit == "0":
                    final_state = "-"
                    final_residue = "NA"
                    final_aa_pp = "NA"
                    aa_source = "not_used_gap"

                else:
                    residue_number += 1
                    final_residue = residue_number

                    if state == "-":
                        final_state = aa_map
                        final_aa_pp = aa_map_pp
                        aa_source = "max_posterior_from_undetermined"
                    else:
                        final_state = state
                        final_aa_pp = record["posteriors"][state]
                        aa_source = "iqtree_state"

                row = {
                    "node": node,
                    "alignment_site": site,
                    "iqtree_state": state,
                    "aa_map": aa_map,
                    "aa_map_pp": f"{aa_map_pp:.6f}",
                    "aa_map_tie": "yes" if len(tied) > 1 else "no",
                    "p_absent": f"{b['p_absent']:.5f}",
                    "p_present": f"{b['p_present']:.5f}",
                    "binary_state": bit,
                    "binary_status": b["status"],
                    "final_state": final_state,
                    "final_residue": final_residue,
                    "final_aa_pp": (
                        final_aa_pp
                        if final_aa_pp == "NA"
                        else f"{final_aa_pp:.6f}"
                    ),
                    "aa_source": aa_source,
                }

                for aa in posterior_labels:
                    col = f"p_{aa}"
                    # Preserve IQ-TREE's original numeric text when available.
                    row[col] = record["raw_posteriors"][col]

                writer.writerow(row)


def write_nogap_posterior_tsv(
    state_records,
    indel_dict,
    binary_posterior,
    posterior_labels,
    outpath,
):
    """
    Write one row per residue in the final gap-stripped ancestral sequence.

    This is the compact protein-oriented table intended for plotting
    residue number vs amino-acid posterior probability and for downstream
    structure/confidence analyses.
    """
    fieldnames = [
        "node",
        "residue",
        "alignment_site",
        "aa",
        "aa_pp",
        "p_present",
        "iqtree_state",
        "aa_source",
    ]

    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for node, records in state_records.items():
            residue_number = 0
            bits = indel_dict[node]

            for record, bit in zip(records, bits):
                if bit == "0":
                    continue

                residue_number += 1
                site = record["site"]
                state = record["state"]

                if state == "-":
                    aa, aa_pp, _ = get_max_posterior_state(
                        record,
                        posterior_labels,
                    )
                    aa_source = "max_posterior_from_undetermined"
                else:
                    if state not in record["posteriors"]:
                        raise ValueError(
                            f"IQ-TREE State '{state}' has no matching posterior "
                            f"column at {node} site {site}."
                        )
                    aa = state
                    aa_pp = record["posteriors"][state]
                    aa_source = "iqtree_state"

                writer.writerow(
                    {
                        "node": node,
                        "residue": residue_number,
                        "alignment_site": site,
                        "aa": aa,
                        "aa_pp": f"{aa_pp:.6f}",
                        "p_present": (
                            f"{binary_posterior[node][site]['p_present']:.5f}"
                        ),
                        "iqtree_state": state,
                        "aa_source": aa_source,
                    }
                )


# ---------------------------
# Resolution-event reporting
# ---------------------------

def report_resolution_events(events):
    if not events:
        print(
            "[INFO] No IQ-TREE State='-' sites required amino-acid "
            "resolution in binary-present positions."
        )
        return

    tie_events = [ev for ev in events if len(ev["tied_aas"]) > 1]

    print(
        f"[INFO] Resolved {len(events)} IQ-TREE State='-' site(s) at "
        "binary-present positions using the highest amino-acid posterior "
        "probability."
    )

    if tie_events:
        print(
            f"[WARN] {len(tie_events)} resolved site(s) had multiple amino "
            "acids tied for the maximum posterior probability.",
            file=sys.stderr,
        )
        print(
            "[WARN] The first tied amino acid in the IQ-TREE .state header "
            "order was retained deterministically.",
            file=sys.stderr,
        )

        for ev in tie_events[:10]:
            tied = ",".join(ev["tied_aas"])
            print(
                f"[WARN]   {ev['node']} site {ev['site']}: "
                f"max_pp={ev['chosen_pp']:.6g}; tied={tied}; "
                f"selected={ev['chosen_aa']}",
                file=sys.stderr,
            )

        if len(tie_events) > 10:
            print(
                f"[WARN]   ... and {len(tie_events) - 10} more tied site(s)",
                file=sys.stderr,
            )


# ---------------------------
# FASTA writer
# ---------------------------

def write_fasta(seqs, outpath, remove_gaps=False):
    with open(outpath, "w") as f:
        for name, seq in seqs.items():
            out_seq = seq.replace("-", "") if remove_gaps else seq
            f.write(f">{name}\n{out_seq}\n")


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine IQ-TREE .state output and RAxML-NG-derived binary "
            "gap-state results into ancestral FASTA files and optional "
            "posterior-probability TSVs."
        )
    )

    parser.add_argument(
        "--state",
        required=True,
        help="IQ-TREE .state file",
    )

    parser.add_argument(
        "--indel",
        required=True,
        help=(
            "Binary gap-state file with IQ-TREE NodeXX labels "
            "(NodeXX followed by a 0/1 state string)"
        ),
    )

    parser.add_argument(
        "--binary_posterior",
        default=None,
        help=(
            "Binary posterior TSV written by raxmlng_binary_asr.py; "
            "required when posterior output TSVs are requested"
        ),
    )

    parser.add_argument(
        "--out_raw_withgap",
        default=None,
        help=(
            "Output raw IQ-TREE ASR FASTA at alignment length "
            "(State column preserved exactly)"
        ),
    )

    parser.add_argument(
        "--out_raw_nogap",
        default=None,
        help=(
            "Output raw IQ-TREE ASR FASTA with '-' removed "
            "(State column otherwise preserved exactly)"
        ),
    )

    parser.add_argument(
        "--out_withgap",
        default="ASR_final_withgap.fasta",
        help="Output indel-aware FASTA with gaps",
    )

    parser.add_argument(
        "--out_nogap",
        default="ASR_final_nogap.fasta",
        help="Output gap-stripped indel-aware FASTA",
    )

    parser.add_argument(
        "--out_site_posterior",
        default=None,
        help=(
            "Output alignment-site posterior TSV including gap-classified sites"
        ),
    )

    parser.add_argument(
        "--out_nogap_posterior",
        default=None,
        help=(
            "Output final gap-stripped residue posterior TSV"
        ),
    )

    args = parser.parse_args()

    posterior_outputs_requested = (
        args.out_site_posterior is not None
        or args.out_nogap_posterior is not None
    )

    if posterior_outputs_requested and args.binary_posterior is None:
        parser.error(
            "--binary_posterior is required when posterior TSV output "
            "is requested."
        )

    try:
        print(f"[INFO] Reading state file: {args.state}")
        raw_state_dict, state_records, posterior_labels = parse_state_file(
            args.state
        )
        print(
            f"[INFO] Loaded {len(raw_state_dict)} sequences from state file"
        )

        print(f"[INFO] Reading binary gap-state file: {args.indel}")
        indel_dict = parse_indel_file(args.indel)
        print(
            f"[INFO] Loaded {len(indel_dict)} binary gap-state profiles"
        )

        validate_node_sets(state_records, indel_dict)

        binary_posterior = None
        if args.binary_posterior:
            print(
                f"[INFO] Reading binary posterior table: "
                f"{args.binary_posterior}"
            )
            binary_posterior = parse_binary_posterior(
                args.binary_posterior
            )
            validate_binary_posterior_against_profiles(
                state_records,
                indel_dict,
                binary_posterior,
            )
            print(
                "[INFO] Binary posterior table matches node/site sets and "
                "binary profiles."
            )

        print("[INFO] Merging state and binary gap-state information...")
        merged, resolution_events = merge_state_and_indel(
            state_records,
            indel_dict,
            posterior_labels,
        )

        print(f"[INFO] Merged {len(merged)} sequences")
        report_resolution_events(resolution_events)

        if len(merged) != len(raw_state_dict):
            raise ValueError(
                f"Merged sequence count ({len(merged)}) does not match "
                f"the number of IQ-TREE ancestral sequences "
                f"({len(raw_state_dict)})."
            )

        if args.out_raw_withgap:
            write_fasta(
                raw_state_dict,
                args.out_raw_withgap,
                remove_gaps=False,
            )

        if args.out_raw_nogap:
            write_fasta(
                raw_state_dict,
                args.out_raw_nogap,
                remove_gaps=True,
            )

        write_fasta(
            merged,
            args.out_withgap,
            remove_gaps=False,
        )

        write_fasta(
            merged,
            args.out_nogap,
            remove_gaps=True,
        )

        if args.out_site_posterior:
            write_site_posterior_tsv(
                state_records,
                indel_dict,
                binary_posterior,
                posterior_labels,
                args.out_site_posterior,
            )
            print(
                f"[DONE] Alignment-site posterior TSV written: "
                f"{args.out_site_posterior}"
            )

        if args.out_nogap_posterior:
            write_nogap_posterior_tsv(
                state_records,
                indel_dict,
                binary_posterior,
                posterior_labels,
                args.out_nogap_posterior,
            )
            print(
                f"[DONE] Final nogap posterior TSV written: "
                f"{args.out_nogap_posterior}"
            )

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

    print("[DONE] Output written:")

    if args.out_raw_withgap:
        print(
            f"           raw with gaps     : "
            f"{args.out_raw_withgap}"
        )

    if args.out_raw_nogap:
        print(
            f"           raw no gaps       : "
            f"{args.out_raw_nogap}"
        )

    print(
        f"           indel with gaps   : {args.out_withgap}"
    )
    print(
        f"           indel no gaps     : {args.out_nogap}"
    )

    if args.out_site_posterior:
        print(
            f"           site posterior    : {args.out_site_posterior}"
        )

    if args.out_nogap_posterior:
        print(
            f"           nogap posterior   : {args.out_nogap_posterior}"
        )


if __name__ == "__main__":
    main()

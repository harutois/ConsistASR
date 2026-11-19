#!/usr/bin/env python3
"""
paml_state_and_indel_to_fasta.py

Combine PAML marginal ancestral amino-acid reconstructions with
RAxML-NG binary presence/absence reconstruction.

The PAML root is retained in raw PAML FASTA and in the complete
alignment-site posterior table, but is excluded from indel-aware FASTA
and final nogap posterior output because RAxML-NG has no one-to-one
binary-ASR counterpart for the degree-2 PAML root.

Outputs
-------
1) Raw PAML ancestral FASTA at alignment length (all PAML nodes)
2) Raw gap-stripped PAML ancestral FASTA (all PAML nodes)
3) Indel-aware FASTA with gaps (mapped non-root nodes only)
4) Indel-aware gap-stripped FASTA (mapped non-root nodes only)
5) Alignment-site posterior TSV (all PAML nodes; root binary fields = NA)
6) Final nogap residue posterior TSV (mapped non-root nodes only)
"""

import argparse
import csv
import re
import sys


BEST_STATE_RE = re.compile(r"([A-Za-z*?-])\((\d*\.?\d+(?:[Ee][+-]?\d+)?)\)")


def parse_sequence_chunk(text, context):
    compact = re.sub(r"\s+", "", text)

    if not compact:
        return ""

    if not re.fullmatch(r"[A-Za-z-]+", compact):
        raise ValueError(
            f"Unexpected character(s) in PAML ancestral sequence "
            f"({context}): '{text.strip()}'"
        )

    return compact.upper()


def parse_all_nodes_from_rst(rst_path):
    """
    Parse ancestral sequences and preserve their node order as written in rst.
    """
    node_chunks = {}
    node_order = []
    current_id = None
    in_section = False
    found_section = False

    node_re = re.compile(
        r"^\s*node\s*#\s*(\d+)\b(.*)$",
        flags=re.IGNORECASE,
    )

    with open(rst_path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            raw = line.rstrip("\r\n")
            stripped = raw.strip()
            lower = stripped.lower()

            if "list of extant and reconstructed sequences" in lower:
                in_section = True
                found_section = True
                current_id = None
                continue

            if not in_section:
                continue

            if (
                lower.startswith("overall")
                or lower.startswith("prob distribution at node")
            ):
                break

            m = node_re.match(raw)
            if m:
                node_id = m.group(1)

                if node_id in node_chunks:
                    raise ValueError(
                        f"Duplicate ancestral node #{node_id} in PAML rst "
                        f"at line {lineno}."
                    )

                current_id = node_id
                node_order.append(node_id)
                node_chunks[current_id] = []

                chunk = parse_sequence_chunk(
                    m.group(2),
                    f"node {node_id}, line {lineno}",
                )

                if chunk:
                    node_chunks[current_id].append(chunk)

                continue

            if current_id is not None:
                if not stripped:
                    current_id = None
                    continue

                chunk = parse_sequence_chunk(
                    raw,
                    f"node {current_id}, line {lineno}",
                )

                if chunk:
                    node_chunks[current_id].append(chunk)

    if not found_section:
        raise ValueError(
            "'List of extant and reconstructed sequences' section not found "
            f"in {rst_path}"
        )

    seqs = {
        node: "".join(node_chunks[node])
        for node in node_order
    }

    if not seqs:
        raise ValueError(f"No ancestral node sequences found in {rst_path}")

    if any(not seq for seq in seqs.values()):
        raise ValueError("Empty PAML ancestral sequence detected.")

    lengths = {len(seq) for seq in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(
            "PAML ancestral sequences have inconsistent lengths: "
            + ", ".join(str(x) for x in sorted(lengths))
        )

    aln_len = next(iter(lengths))
    print(f"[INFO] PAML ancestral alignment length: {aln_len}")

    return seqs, node_order, aln_len


def parse_best_state_pp_from_rst(rst_path, node_order, alignment_length):
    """
    Parse PAML block:
        Prob of best state at each node, listed by site

    Each site row contains one X(PP) entry per ancestral node.
    The entry order is matched to the node order in the reconstructed
    sequence section of the same rst file.
    """
    with open(rst_path, encoding="utf-8") as f:
        lines = f.read().splitlines()

    marker = None
    for i, line in enumerate(lines):
        if "prob of best state at each node, listed by site" in line.lower():
            marker = i
            break

    if marker is None:
        raise ValueError(
            "'Prob of best state at each node, listed by site' block "
            "not found in PAML rst."
        )

    header_seen = False
    records = {node: {} for node in node_order}
    site_re = re.compile(r"^\s*(\d+)\s+\d+\s+.*?:\s*(.*)$")

    for line in lines[marker + 1:]:
        lower = line.lower()

        if not header_seen:
            if "site" in lower and "freq" in lower and "data" in lower:
                header_seen = True
            continue

        m = site_re.match(line)
        if not m:
            # Once site rows have started, stop at the next major section.
            if records[node_order[0]]:
                stripped = line.strip()
                if stripped and (
                    stripped.startswith("(")
                    or "list of extant and reconstructed sequences" in lower
                    or lower.startswith("prob distribution at node")
                ):
                    break
            continue

        site = int(m.group(1))
        entries = BEST_STATE_RE.findall(m.group(2))

        if len(entries) != len(node_order):
            raise ValueError(
                f"PAML best-state PP row at site {site} contains "
                f"{len(entries)} ancestral entries; expected "
                f"{len(node_order)}."
            )

        for node, (state, pp_text) in zip(node_order, entries):
            pp = float(pp_text)

            if pp < 0.0 or pp > 1.0:
                raise ValueError(
                    f"PAML posterior probability out of range at "
                    f"node {node} site {site}: {pp}"
                )

            records[node][site] = {
                "state": state.upper(),
                "pp": pp,
                "pp_text": pp_text,
            }

    expected_sites = set(range(1, alignment_length + 1))

    for node in node_order:
        found = set(records[node])

        if found != expected_sites:
            raise ValueError(
                f"PAML best-state PP sites for node {node} do not match "
                f"1-{alignment_length}; missing="
                f"{sorted(expected_sites - found)[:10]}"
            )

    return records


def validate_paml_state_vs_pp(seqs, pp_records):
    for node, seq in seqs.items():
        for site, aa in enumerate(seq, start=1):
            pp_state = pp_records[node][site]["state"]

            if aa != pp_state:
                raise ValueError(
                    f"PAML reconstructed sequence / best-state PP mismatch at "
                    f"node {node} site {site}: sequence={aa}, PP_state={pp_state}"
                )


def parse_indel_file(indel_file):
    indels = {}

    with open(indel_file, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip():
                continue

            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Malformed binary gap-state line {lineno}: "
                    f"{line.strip()}"
                )

            node, bits = parts

            if node in indels:
                raise ValueError(
                    f"Duplicate PAML node in binary gap-state file: {node}"
                )

            if not re.fullmatch(r"[01]+", bits):
                raise ValueError(
                    f"Invalid binary profile for PAML node {node}; "
                    "only 0/1 allowed."
                )

            indels[node] = list(bits)

    if not indels:
        raise ValueError("No mapped binary gap-state profiles found.")

    return indels


def parse_binary_posterior(path):
    required = {
        "paml_node",
        "raxmlng_node",
        "site",
        "raxmlng_state",
        "p_absent",
        "p_present",
        "map_state",
        "merge_state",
        "status",
    }

    data = {}

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("Mapped binary posterior TSV has no header.")

        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Mapped binary posterior TSV missing column(s): "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            node = row["paml_node"]
            site = int(row["site"])
            key = (node, site)

            if key in data:
                raise ValueError(
                    f"Duplicate mapped binary posterior row: "
                    f"node {node} site {site}"
                )

            p_absent = float(row["p_absent"])
            p_present = float(row["p_present"])

            if not (
                0.0 <= p_absent <= 1.0
                and 0.0 <= p_present <= 1.0
            ):
                raise ValueError(
                    f"Mapped binary posterior out of range at "
                    f"node {node} site {site}."
                )

            data[key] = {
                **row,
                "p_absent_float": p_absent,
                "p_present_float": p_present,
            }

    if not data:
        raise ValueError("No rows in mapped binary posterior TSV.")

    return data


def validate_node_sets_and_lengths(
    seqs,
    indels,
    binary_pp,
    root_node,
    alignment_length,
):
    seq_nodes = set(seqs)
    indel_nodes = set(indels)

    expected_nonroot = seq_nodes - {root_node}

    if root_node not in seq_nodes:
        raise ValueError(
            f"Declared PAML root node {root_node} is absent from rst sequences."
        )

    if indel_nodes != expected_nonroot:
        raise ValueError(
            "Mapped binary node set must equal all PAML ancestral nodes "
            "except the root; "
            f"missing={sorted(expected_nonroot - indel_nodes, key=int)}, "
            f"extra={sorted(indel_nodes - expected_nonroot)}"
        )

    pp_nodes = {node for node, _ in binary_pp}
    if pp_nodes != expected_nonroot:
        raise ValueError(
            "Mapped binary posterior node set does not equal PAML non-root "
            "node set."
        )

    for node in expected_nonroot:
        if len(indels[node]) != alignment_length:
            raise ValueError(
                f"Binary profile length mismatch for PAML node {node}: "
                f"{len(indels[node])} vs {alignment_length}"
            )

        for site in range(1, alignment_length + 1):
            row = binary_pp.get((node, site))
            if row is None:
                raise ValueError(
                    f"Missing binary posterior row for node {node} site {site}."
                )

            if row["merge_state"] != indels[node][site - 1]:
                raise ValueError(
                    f"Binary profile/posterior mismatch at node {node} "
                    f"site {site}."
                )


def merge_nonroot_sequences(seqs, indels, root_node):
    merged = {}

    for node, seq in seqs.items():
        if node == root_node:
            continue

        bits = indels[node]
        merged[node] = "".join(
            aa if bit == "1" else "-"
            for aa, bit in zip(seq, bits)
        )

    return merged


def write_fasta(seqs, outpath, remove_gaps=False):
    with open(outpath, "w", encoding="utf-8") as f:
        for node, seq in seqs.items():
            out_seq = seq.replace("-", "") if remove_gaps else seq
            f.write(f">{node}\n{out_seq}\n")


def write_site_posterior(
    seqs,
    node_order,
    paml_pp,
    indels,
    binary_pp,
    root_node,
    outpath,
):
    fields = [
        "node",
        "alignment_site",
        "paml_state",
        "aa_pp",
        "p_absent",
        "p_present",
        "binary_state",
        "binary_status",
        "final_state",
        "final_residue",
        "integration_status",
    ]

    with open(outpath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for node in node_order:
            residue_number = 0
            seq = seqs[node]

            for site, aa in enumerate(seq, start=1):
                pp = paml_pp[node][site]["pp_text"]

                if node == root_node:
                    writer.writerow(
                        {
                            "node": node,
                            "alignment_site": site,
                            "paml_state": aa,
                            "aa_pp": pp,
                            "p_absent": "NA",
                            "p_present": "NA",
                            "binary_state": "NA",
                            "binary_status": "NA",
                            "final_state": "NA",
                            "final_residue": "NA",
                            "integration_status": "paml_root_not_processed",
                        }
                    )
                    continue

                bit = indels[node][site - 1]
                b = binary_pp[(node, site)]

                if bit == "1":
                    residue_number += 1
                    final_state = aa
                    final_residue = residue_number
                    integration_status = "present"
                else:
                    final_state = "-"
                    final_residue = "NA"
                    integration_status = "absent"

                writer.writerow(
                    {
                        "node": node,
                        "alignment_site": site,
                        "paml_state": aa,
                        "aa_pp": pp,
                        "p_absent": b["p_absent"],
                        "p_present": b["p_present"],
                        "binary_state": bit,
                        "binary_status": b["status"],
                        "final_state": final_state,
                        "final_residue": final_residue,
                        "integration_status": integration_status,
                    }
                )


def write_nogap_posterior(
    seqs,
    node_order,
    paml_pp,
    indels,
    binary_pp,
    root_node,
    outpath,
):
    fields = [
        "node",
        "residue",
        "alignment_site",
        "aa",
        "aa_pp",
        "p_present",
    ]

    with open(outpath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for node in node_order:
            if node == root_node:
                continue

            residue_number = 0

            for site, aa in enumerate(seqs[node], start=1):
                if indels[node][site - 1] == "0":
                    continue

                residue_number += 1

                writer.writerow(
                    {
                        "node": node,
                        "residue": residue_number,
                        "alignment_site": site,
                        "aa": aa,
                        "aa_pp": paml_pp[node][site]["pp_text"],
                        "p_present": binary_pp[(node, site)]["p_present"],
                    }
                )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rst", required=True)
    parser.add_argument("--indel", required=True)
    parser.add_argument("--binary-posterior", required=True)
    parser.add_argument("--root-node", required=True)
    parser.add_argument("--out_raw_withgap", default=None)
    parser.add_argument("--out_raw_nogap", default=None)
    parser.add_argument(
        "--out_withgap",
        default="ASR_PAML_indel_withgap.fasta",
    )
    parser.add_argument(
        "--out_nogap",
        default="ASR_PAML_indel_nogap.fasta",
    )
    parser.add_argument("--out-site-posterior", default=None)
    parser.add_argument("--out-nogap-posterior", default=None)
    args = parser.parse_args()

    try:
        seqs, node_order, aln_len = parse_all_nodes_from_rst(args.rst)
        print(f"[INFO] Loaded {len(seqs)} ancestral nodes from rst")

        paml_pp = parse_best_state_pp_from_rst(
            args.rst,
            node_order,
            aln_len,
        )
        validate_paml_state_vs_pp(seqs, paml_pp)
        print(
            "[INFO] PAML best-state posterior probabilities match "
            "reconstructed sequences."
        )

        indels = parse_indel_file(args.indel)
        print(
            f"[INFO] Loaded {len(indels)} mapped non-root binary profiles"
        )

        binary_pp = parse_binary_posterior(args.binary_posterior)

        validate_node_sets_and_lengths(
            seqs,
            indels,
            binary_pp,
            args.root_node,
            aln_len,
        )

        print(
            f"[INFO] PAML root node {args.root_node} retained only in raw "
            "PAML outputs and site-level posterior table."
        )

        merged = merge_nonroot_sequences(
            seqs,
            indels,
            args.root_node,
        )
        print(
            f"[INFO] Merged {len(merged)} non-root ancestral sequences"
        )

        if args.out_raw_withgap:
            write_fasta(
                seqs,
                args.out_raw_withgap,
                remove_gaps=False,
            )

        if args.out_raw_nogap:
            write_fasta(
                seqs,
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
            write_site_posterior(
                seqs,
                node_order,
                paml_pp,
                indels,
                binary_pp,
                args.root_node,
                args.out_site_posterior,
            )
            print(
                f"[DONE] Alignment-site posterior TSV written: "
                f"{args.out_site_posterior}"
            )

        if args.out_nogap_posterior:
            write_nogap_posterior(
                seqs,
                node_order,
                paml_pp,
                indels,
                binary_pp,
                args.root_node,
                args.out_nogap_posterior,
            )
            print(
                f"[DONE] Final nogap posterior TSV written: "
                f"{args.out_nogap_posterior}"
            )

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print("[DONE] Output written:")

    if args.out_raw_withgap:
        print(f"           raw with gaps     : {args.out_raw_withgap}")
    if args.out_raw_nogap:
        print(f"           raw no gaps       : {args.out_raw_nogap}")

    print(f"           indel with gaps   : {args.out_withgap}")
    print(f"           indel no gaps     : {args.out_nogap}")

    if args.out_site_posterior:
        print(f"           site posterior    : {args.out_site_posterior}")
    if args.out_nogap_posterior:
        print(f"           nogap posterior   : {args.out_nogap_posterior}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3

import argparse
import math
import sys

import gemmi


# 3-letter to 1-letter amino acid code
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


def read_cif_sequence_and_atoms(cif_path, chain_id=None):
    """
    Read an AlphaFold mmCIF file and return:
      - structure sequence in 1-letter codes,
      - a list of atom lists, one per residue,
      - the gemmi.Structure object,
      - the actual chain ID used.

    If --chain is not specified, exactly one polymer chain is expected.
    """

    try:
        st = gemmi.read_structure(cif_path)
    except Exception as exc:
        sys.exit(f"[ERROR] Failed to read mmCIF file {cif_path}: {exc}")

    if len(st) == 0:
        sys.exit(f"[ERROR] No model found in mmCIF file: {cif_path}")

    mdl = st[0]
    chains = list(mdl)

    def polymer_residues(ch):
        return [
            res for res in ch
            if getattr(res, "entity_type", None) == gemmi.EntityType.Polymer
        ]

    # -------------------------
    # Chain selection
    # -------------------------
    if chain_id is not None:
        matches = [
            ch for ch in chains
            if ch.name.strip() == chain_id.strip()
        ]

        if len(matches) == 0:
            sys.exit(
                f"[ERROR] Chain '{chain_id}' not found in structure."
            )

        if len(matches) > 1:
            sys.exit(
                f"[ERROR] Multiple chains named '{chain_id}' found "
                f"in structure."
            )

        ch = matches[0]

        if not polymer_residues(ch):
            sys.exit(
                f"[ERROR] Chain '{chain_id}' is not a polymer chain."
            )

    else:
        polys = [
            ch for ch in chains
            if polymer_residues(ch)
        ]

        if not polys:
            sys.exit("[ERROR] No polymer chain found in structure.")

        if len(polys) > 1:
            desc = ", ".join(
                f"{ch.name!r} ({len(polymer_residues(ch))} residues)"
                for ch in polys
            )

            sys.exit(
                "[ERROR] Multiple polymer chains found in structure: "
                f"{desc}. Specify the intended chain with --chain."
            )

        ch = polys[0]
        chain_id = ch.name

    # -------------------------
    # Build structure sequence
    # -------------------------
    seq = []
    res_atoms = []

    for res in polymer_residues(ch):
        aa3 = res.name.strip().upper()
        aa1 = AA3_TO_1.get(aa3, "X")

        atoms = list(res)

        if not atoms:
            sys.exit(
                f"[ERROR] Residue {res.name} {res.seqid} "
                f"in chain {chain_id!r} contains no atoms."
            )

        if not any(atom.name.strip() == "CA" for atom in atoms):
            sys.exit(
                f"[ERROR] Residue {res.name} {res.seqid} "
                f"in chain {chain_id!r} has no CA atom."
            )

        seq.append(aa1)
        res_atoms.append(atoms)

    if not seq:
        sys.exit(
            f"[ERROR] No polymer residues found in chain {chain_id!r}."
        )

    return "".join(seq), res_atoms, st, chain_id


def read_single_fasta(path):
    """
    Read exactly one FASTA entry.

    Multi-line sequences are supported.
    Multiple FASTA entries are rejected to avoid ambiguous mapping.
    """

    entries = []
    header = None
    seq_parts = []

    def store_current():
        if header is not None:
            entries.append((header, "".join(seq_parts)))

    with open(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                store_current()
                header = line[1:]
                seq_parts = []
            else:
                if header is None:
                    sys.exit(
                        f"[ERROR] Sequence data found before the first "
                        f"FASTA header at line {lineno} in {path}."
                    )

                seq_parts.append("".join(line.split()))

    store_current()

    if not entries:
        sys.exit(f"[ERROR] No FASTA entry found in {path}.")

    if len(entries) != 1:
        sys.exit(
            f"[ERROR] Expected exactly one FASTA entry in {path}, "
            f"but found {len(entries)}."
        )

    header, sequence = entries[0]

    if not sequence:
        sys.exit(
            f"[ERROR] Empty FASTA sequence for >{header} in {path}."
        )

    sequence = sequence.upper()

    bad = sorted(
        set(c for c in sequence if not (c.isalpha() or c == "-"))
    )

    if bad:
        sys.exit(
            f"[ERROR] Unsupported character(s) in ASR FASTA: "
            f"{', '.join(repr(x) for x in bad)}"
        )

    return sequence


def read_pp_values(path):
    """Read per-site PP values, one numeric value per non-empty line."""

    vals = []

    with open(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            s = raw.strip()

            if not s:
                continue

            try:
                value = float(s)
            except ValueError:
                sys.exit(
                    f"[ERROR] Non-numeric PP value at line {lineno} "
                    f"in {path}: {s!r}"
                )

            if not math.isfinite(value):
                sys.exit(
                    f"[ERROR] Non-finite PP value at line {lineno} "
                    f"in {path}: {s!r}"
                )

            vals.append(value)

    if not vals:
        sys.exit(f"[ERROR] No PP values found in {path}.")

    return vals


def make_ungapped_to_aln_map(aln_seq):
    """
    From an aligned sequence containing '-' gaps, return:
      - ungapped sequence,
      - mapping from ungapped residue index to alignment column.
    """

    map_list = []
    ungapped = []

    for i, char in enumerate(aln_seq):
        if char != "-":
            map_list.append(i)
            ungapped.append(char)

    return "".join(ungapped), map_list


def report_seq_mismatch(name_a, seq_a, name_b, seq_b, max_show=10):
    """Return a compact sequence-mismatch description."""

    diffs = []

    for i, (aa, bb) in enumerate(
        zip(seq_a, seq_b),
        start=1,
    ):
        if aa != bb:
            diffs.append((i, aa, bb))

            if len(diffs) >= max_show:
                break

    parts = []

    if diffs:
        parts.append(
            "first difference(s): "
            + ", ".join(
                f"{pos}:{aa}->{bb}"
                for pos, aa, bb in diffs
            )
        )

    if len(seq_a) != len(seq_b):
        parts.append(
            f"lengths: {name_a}={len(seq_a)}, "
            f"{name_b}={len(seq_b)}"
        )

    if not parts:
        parts.append("sequence mismatch detected")

    return "; ".join(parts)


def scale_pp_to_100(pp_vals):
    """
    Scale PP values to 0-100.

    If max(PP) <= 1.00001, values are interpreted as probabilities
    on the 0-1 scale. Otherwise they are interpreted as already being
    on the 0-100 scale.
    """

    if any(v < 0.0 for v in pp_vals):
        sys.exit("[ERROR] Negative PP value detected.")

    mx = max(pp_vals)

    if mx <= 1.00001:
        scaled = [v * 100.0 for v in pp_vals]
        info = "0-1 -> x100"

    elif mx <= 100.0:
        scaled = list(pp_vals)
        info = "already 0-100"

    else:
        sys.exit(
            f"[ERROR] PP value exceeds expected range: max={mx}"
        )

    for value in scaled:
        if not (0.0 <= value <= 100.001):
            sys.exit(
                f"[ERROR] Scaled PP value outside 0-100 range: "
                f"{value}"
            )

    return scaled, info


def per_residue_plddt(res_atoms):
    """
    Extract one pLDDT value per residue from the CA atom.

    AlphaFold stores residue-level pLDDT in the B-factor field.
    """

    values = []

    for i, atoms in enumerate(res_atoms, start=1):
        ca_atoms = [
            atom for atom in atoms
            if atom.name.strip() == "CA"
        ]

        if len(ca_atoms) != 1:
            sys.exit(
                f"[ERROR] Expected exactly one CA atom for "
                f"structure residue {i}, found {len(ca_atoms)}."
            )

        value = ca_atoms[0].b_iso

        if not math.isfinite(value):
            sys.exit(
                f"[ERROR] Non-finite pLDDT value at "
                f"structure residue {i}."
            )

        if not (0.0 <= value <= 100.0):
            sys.exit(
                f"[ERROR] pLDDT outside 0-100 range at "
                f"structure residue {i}: {value}"
            )

        values.append(value)

    return values


def write_pdb_with_b(st, res_atoms, bvals, out_path):
    """
    Write a PDB file where all atoms of residue i receive
    the same B-factor bvals[i].
    """

    if len(res_atoms) != len(bvals):
        sys.exit(
            f"[ERROR] Number of residues ({len(res_atoms)}) "
            f"does not match number of B-factor values "
            f"({len(bvals)})."
        )

    for i, (atoms, bval) in enumerate(
        zip(res_atoms, bvals),
        start=1,
    ):
        if not math.isfinite(bval):
            sys.exit(
                f"[ERROR] Non-finite B-factor value at residue {i}."
            )

        for atom in atoms:
            atom.b_iso = bval

    try:
        with open(out_path, "w") as fh:
            fh.write(st.make_pdb_string())
    except Exception as exc:
        sys.exit(
            f"[ERROR] Failed to write PDB file {out_path}: {exc}"
        )

    print(f"[OK] Wrote: {out_path}")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Map ASR site-wise PP and/or AlphaFold pLDDT "
            "to PDB B-factors while preserving alignment-gap "
            "correspondence."
        )
    )

    ap.add_argument(
        "--cif",
        required=True,
        help="AlphaFold mmCIF file",
    )

    ap.add_argument(
        "--asr_fasta",
        required=True,
        help=(
            "Single aligned ancestral FASTA entry "
            "containing '-' gaps"
        ),
    )

    ap.add_argument(
        "--pp_file",
        required=True,
        help=(
            "Per-site PP values file "
            "(one value per alignment column)"
        ),
    )

    ap.add_argument(
        "--chain",
        default=None,
        help=(
            "Chain ID in structure. If omitted, the structure "
            "must contain exactly one polymer chain."
        ),
    )

    ap.add_argument(
        "--mode",
        choices=["plddt", "pp", "diff", "prod"],
        required=True,
        help=(
            "Metric written to B-factor: "
            "plddt = AlphaFold pLDDT; "
            "pp = ASR posterior probability; "
            "diff = scaled PP-pLDDT difference "
            "((PP-pLDDT+100)/2); "
            "prod = PPxpLDDT/100"
        ),
    )

    ap.add_argument(
        "--out",
        required=True,
        help="Output PDB path",
    )

    args = ap.parse_args()

    # --------------------------------------------------------
    # 1) Read structure
    # --------------------------------------------------------
    struct_seq, res_atoms, st, chain_id = (
        read_cif_sequence_and_atoms(
            args.cif,
            args.chain,
        )
    )

    print(
        f"[INFO] Structure chain={chain_id}, "
        f"residues={len(res_atoms)}"
    )

    # --------------------------------------------------------
    # 2) Read ancestral sequence and PP
    # --------------------------------------------------------
    asr_aln = read_single_fasta(args.asr_fasta)
    pp_vals = read_pp_values(args.pp_file)

    if len(pp_vals) != len(asr_aln):
        sys.exit(
            f"[ERROR] PP length ({len(pp_vals)}) != "
            f"alignment length ({len(asr_aln)}). "
            "The PP file and aligned ancestral FASTA must "
            "use the same alignment-coordinate system."
        )

    # --------------------------------------------------------
    # 3) Build alignment -> residue mapping
    # --------------------------------------------------------
    asr_ungap, map_u2aln = make_ungapped_to_aln_map(
        asr_aln
    )

    print(
        f"[INFO] Alignment length={len(asr_aln)}, "
        f"ungapped length={len(asr_ungap)}"
    )

    # --------------------------------------------------------
    # 4) Require exact ancestral/structure sequence identity
    # --------------------------------------------------------
    if asr_ungap != struct_seq:
        detail = report_seq_mismatch(
            "ASR(ungapped)",
            asr_ungap,
            "Structure",
            struct_seq,
        )

        sys.exit(
            "[ERROR] Ungapped ancestral sequence does not "
            f"exactly match the structure sequence; {detail}."
        )

    print(
        "[INFO] Structure sequence exactly matches "
        "the ungapped ancestral sequence."
    )

    # --------------------------------------------------------
    # 5) Map PP to retained ancestral residues
    # --------------------------------------------------------
    pp_scaled, scale_info = scale_pp_to_100(pp_vals)

    print(f"[INFO] PP scale: {scale_info}")

    per_res_pp = [
        pp_scaled[aln_col]
        for aln_col in map_u2aln
    ]

    if len(per_res_pp) != len(res_atoms):
        sys.exit(
            f"[ERROR] Residue-count mismatch after gap mapping: "
            f"PP residues={len(per_res_pp)}, "
            f"structure residues={len(res_atoms)}."
        )

    # --------------------------------------------------------
    # 6) Extract per-residue AlphaFold pLDDT
    # --------------------------------------------------------
    per_res_plddt = per_residue_plddt(res_atoms)

    if len(per_res_plddt) != len(per_res_pp):
        sys.exit(
            f"[ERROR] PP/pLDDT residue-count mismatch: "
            f"PP={len(per_res_pp)}, "
            f"pLDDT={len(per_res_plddt)}."
        )

    # --------------------------------------------------------
    # 7) Select metric
    # --------------------------------------------------------
    if args.mode == "pp":
        bvals = per_res_pp

    elif args.mode == "plddt":
        bvals = per_res_plddt

    elif args.mode == "diff":
        # PP - pLDDT ranges from -100 to +100.
        # Store a scaled 0-100 value:
        #
        #   (PP - pLDDT + 100) / 2
        #
        bvals = [
            (pp - pl + 100.0) / 2.0
            for pp, pl in zip(
                per_res_pp,
                per_res_plddt,
            )
        ]

    else:  # prod
        # PP x pLDDT / 100 -> 0-100
        bvals = [
            (pp * pl) / 100.0
            for pp, pl in zip(
                per_res_pp,
                per_res_plddt,
            )
        ]

    if len(bvals) != len(res_atoms):
        sys.exit(
            "[ERROR] Internal metric/residue length mismatch."
        )

    # --------------------------------------------------------
    # 8) Write output PDB
    # --------------------------------------------------------
    write_pdb_with_b(
        st,
        res_atoms,
        bvals,
        args.out,
    )


if __name__ == "__main__":
    main()
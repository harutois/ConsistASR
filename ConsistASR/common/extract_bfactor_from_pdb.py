#!/usr/bin/env python3
# ============================================================
# extract_bfactor_from_pdb.py
#
# Extract B-factor values from PDB ATOM records.
#
# Default behavior:
#   - extracts CA atoms
#   - writes a TSV table with residue metadata and B-factor values
#
# Example:
#   python3 extract_bfactor_from_pdb.py \
#     --in Node10_pp_bfactor.pdb \
#     --out Node10_pp_bfactor.tsv
#
# Multiple PDB files:
#   python3 extract_bfactor_from_pdb.py \
#     --in Node10_pp_bfactor.pdb Node10_plddt_bfactor.pdb \
#     --out Node10_bfactor_values.tsv
#
# Values only, for simple downstream plotting:
#   python3 extract_bfactor_from_pdb.py \
#     --in Node10_pp_bfactor.pdb \
#     --values-only > Node10_pp_values.txt
# ============================================================

import argparse
import math
import os
import sys


def parse_pdb_bfactors(pdb_path, atom_name="CA", chain_id=None):
    """
    Extract B-factor values from PDB ATOM records.

    Parameters
    ----------
    pdb_path : str
        Input PDB file.
    atom_name : str
        Atom name to extract, e.g. "CA".
    chain_id : str or None
        Optional chain ID filter.

    Returns
    -------
    list[dict]
        One record per matching atom.
    """

    records = []
    current_model = ""

    try:
        fh = open(pdb_path, encoding="utf-8")
    except OSError as exc:
        sys.exit(
            f"[ERROR] Failed to open PDB file {pdb_path}: {exc}"
        )

    with fh:
        for lineno, line in enumerate(fh, start=1):

            if line.startswith("MODEL"):
                parts = line.split()
                current_model = parts[1] if len(parts) > 1 else ""
                continue

            if not line.startswith("ATOM"):
                continue

            # PDB B-factor field occupies columns 61-66.
            # A valid ATOM record must therefore contain at least 66 chars.
            if len(line.rstrip("\n\r")) < 66:
                sys.exit(
                    f"[ERROR] Malformed ATOM record in {pdb_path} "
                    f"at line {lineno}: line is too short."
                )

            atom = line[12:16].strip()

            if atom != atom_name:
                continue

            chain = line[21].strip()

            if chain_id is not None and chain != chain_id:
                continue

            raw_bfactor = line[60:66]

            try:
                bfactor = float(raw_bfactor)
            except ValueError:
                sys.exit(
                    f"[ERROR] Invalid B-factor in {pdb_path} "
                    f"at line {lineno}: {raw_bfactor!r}"
                )

            if not math.isfinite(bfactor):
                sys.exit(
                    f"[ERROR] Non-finite B-factor in {pdb_path} "
                    f"at line {lineno}: {bfactor}"
                )

            records.append(
                {
                    "file": os.path.basename(pdb_path),
                    "model": current_model,
                    "chain": chain,
                    "resseq": line[22:26].strip(),
                    "icode": line[26].strip(),
                    "resname": line[17:20].strip(),
                    "atom": atom,
                    "bfactor": bfactor,
                }
            )

    return records


def write_tsv(records, out_handle, values_only=False, header=True):
    """
    Write extracted records as TSV.

    If values_only=True, write only one B-factor value per line.
    """

    if values_only:
        for rec in records:
            out_handle.write(f"{rec['bfactor']:.2f}\n")
        return

    columns = [
        "file",
        "model",
        "chain",
        "resseq",
        "icode",
        "resname",
        "atom",
        "bfactor",
    ]

    if header:
        out_handle.write("\t".join(columns) + "\n")

    for rec in records:
        row = []

        for col in columns:
            if col == "bfactor":
                row.append(f"{rec[col]:.2f}")
            else:
                row.append(str(rec[col]))

        out_handle.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract B-factor values from one or more PDB files."
        )
    )

    parser.add_argument(
        "--in",
        dest="input_pdbs",
        required=True,
        nargs="+",
        help="Input PDB file(s)",
    )

    parser.add_argument(
        "--out",
        default=None,
        help="Output TSV file. If omitted, write to stdout.",
    )

    parser.add_argument(
        "--atom",
        default="CA",
        help="Atom name to extract. Default: CA",
    )

    parser.add_argument(
        "--chain",
        default=None,
        help="Optional chain ID filter, e.g. A",
    )

    parser.add_argument(
        "--values-only",
        action="store_true",
        help="Write only B-factor values, one per line.",
    )

    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write a header line in TSV output.",
    )

    args = parser.parse_args()

    all_records = []

    for pdb_path in args.input_pdbs:

        if not os.path.isfile(pdb_path):
            print(
                f"[ERROR] File not found: {pdb_path}",
                file=sys.stderr,
            )
            sys.exit(1)

        records = parse_pdb_bfactors(
            pdb_path,
            atom_name=args.atom,
            chain_id=args.chain,
        )

        if not records:
            print(
                f"[ERROR] No matching ATOM records found in {pdb_path} "
                f"(atom={args.atom}, chain={args.chain or 'any'})",
                file=sys.stderr,
            )
            sys.exit(1)

        all_records.extend(records)

    if not all_records:
        print(
            "[ERROR] No B-factor records were extracted.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.out:
        try:
            with open(
                args.out,
                "w",
                encoding="utf-8",
            ) as out:
                write_tsv(
                    all_records,
                    out,
                    values_only=args.values_only,
                    header=not args.no_header,
                )
        except OSError as exc:
            print(
                f"[ERROR] Failed to write output file "
                f"{args.out}: {exc}",
                file=sys.stderr,
            )
            sys.exit(1)

    else:
        write_tsv(
            all_records,
            sys.stdout,
            values_only=args.values_only,
            header=not args.no_header,
        )


if __name__ == "__main__":
    main()
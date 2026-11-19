#!/usr/bin/env python3
"""
run_paml_asr.py

Minimal helper to:
  - prepare a PAML-compatible tree file
  - generate a codeml .ctl file for amino-acid ASR
  - run codeml
  - keep all PAML outputs in a specified directory

This is primarily an example helper for reproducing the toy 7TM
rhodopsin PAML ASR used in the ConsistASR examples. It also permits
an alternative empirical amino-acid rate file and PAML AA model 2
(Empirical) or 3 (Empirical+F).
"""

import argparse
import os
import subprocess
import sys
import textwrap
from io import StringIO
from pathlib import Path
from shutil import which

from Bio import Phylo


SCRIPT_DIR = Path(__file__).resolve().parent
BUNDLED_LG = SCRIPT_DIR / "ASR" / "lg.dat"


def count_fasta_sequences(fasta_path):
    """Count sequences in a FASTA file."""
    nseq = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                nseq += 1

    if nseq == 0:
        sys.exit(f"[ERROR] No FASTA sequences found in: {fasta_path}")

    return nseq


def prepare_paml_tree(input_tree, output_tree, expected_taxa):
    """
    Prepare a PAML-compatible tree file.

    PAML accepts a header specifying:
        <number_of_taxa> <number_of_trees>

    Example:
        18 1
        ((A:0.1,B:0.1):0.2,C:0.3);
    """

    with open(input_tree, "r") as f:
        tree_text = f.read().strip()

    if not tree_text:
        sys.exit(f"[ERROR] Tree file is empty: {input_tree}")

    if not tree_text.endswith(";"):
        sys.exit(
            f"[ERROR] Newick tree does not end with ';': {input_tree}"
        )

    # Validate the Newick tree and count terminal taxa.
    try:
        tree = Phylo.read(StringIO(tree_text), "newick")
    except Exception as e:
        sys.exit(f"[ERROR] Failed to parse Newick tree: {e}")

    tree_taxa = len(tree.get_terminals())

    if tree_taxa != expected_taxa:
        sys.exit(
            "[ERROR] Number of taxa differs between alignment and tree:\n"
            f"        alignment: {expected_taxa}\n"
            f"        tree     : {tree_taxa}"
        )

    with open(output_tree, "w") as f:
        f.write(f"{expected_taxa} 1\n")
        f.write(tree_text)
        f.write("\n")

    print(f"[INFO] Alignment sequences: {expected_taxa}")
    print(f"[INFO] Tree terminal taxa : {tree_taxa}")
    print(f"[INFO] Wrote PAML-compatible tree: {output_tree}")


def main():
    ap = argparse.ArgumentParser(
        description="Generate a PAML codeml ctl file and run codeml for AA ASR."
    )

    ap.add_argument(
        "--msa", "--aln",
        dest="msa",
        required=True,
        help=(
            "Input amino-acid alignment file used as seqfile in codeml ctl. "
            "--aln is accepted as a backward-compatible alias."
        ),
    )

    ap.add_argument(
        "--tree",
        required=True,
        help="Input tree file in Newick format.",
    )

    ap.add_argument(
        "--outdir",
        default="paml_asr",
        help="Output directory for ctl / codeml outputs (default: paml_asr).",
    )

    ap.add_argument(
        "--prefix",
        default="PAML_LG_FG",
        help="Prefix used for ctl / result files (default: PAML_LG_FG).",
    )

    ap.add_argument(
        "--model",
        type=int,
        choices=(2, 3),
        default=3,
        help=(
            "PAML amino-acid model: 2=Empirical, 3=Empirical+F "
            "(default: 3)."
        ),
    )

    ap.add_argument(
        "--aa-ratefile",
        default=str(BUNDLED_LG),
        help=(
            "Empirical amino-acid rate file used by codeml. "
            "By default, the bundled toy-example ASR/lg.dat is used. "
            "An alternative PAML-compatible rate file may be supplied."
        ),
    )

    args = ap.parse_args()

    # Resolve input paths before changing the working directory for codeml.
    args.msa = os.path.abspath(args.msa)
    args.tree = os.path.abspath(args.tree)
    args.aa_ratefile = os.path.abspath(args.aa_ratefile)
    args.outdir = os.path.abspath(args.outdir)

    # Basic input checks.
    for path in [args.msa, args.tree, args.aa_ratefile]:
        if not os.path.isfile(path):
            sys.exit(f"[ERROR] File not found: {path}")

    if not which("codeml"):
        sys.exit("[ERROR] codeml not found in PATH.")

    os.makedirs(args.outdir, exist_ok=True)

    ctl_path = os.path.join(
        args.outdir,
        f"{args.prefix}.ctl",
    )

    out_path = os.path.join(
        args.outdir,
        f"{args.prefix}.out",
    )

    paml_tree_path = os.path.join(
        args.outdir,
        f"{args.prefix}.paml_tree.nwk",
    )

    # Count alignment sequences.
    nseq = count_fasta_sequences(args.msa)

    # Prepare a PAML-compatible tree with "<n_taxa> 1" header.
    prepare_paml_tree(
        args.tree,
        paml_tree_path,
        nseq,
    )

    # codeml is executed with cwd=args.outdir.
    # Therefore, paths in the ctl file are expressed relative to outdir.
    msa_rel = os.path.relpath(
        os.path.abspath(args.msa),
        args.outdir,
    )

    aa_rate_rel = os.path.relpath(
        os.path.abspath(args.aa_ratefile),
        args.outdir,
    )

    # The PAML-formatted tree is already inside outdir.
    tree_rel = os.path.basename(paml_tree_path)

    ctl_text = f"""
    seqfile = {msa_rel}
    treefile = {tree_rel}
    outfile = {os.path.basename(out_path)}

    noisy = 9
    verbose = 1
    runmode = 0

    seqtype = 2

    model = {args.model}
    aaRatefile = {aa_rate_rel}

    fix_alpha = 0
    alpha = 0.5
    ncatG = 4

    RateAncestor = 1
    cleandata = 0
    fix_blength = 0
    """

    with open(ctl_path, "w") as f:
        f.write(textwrap.dedent(ctl_text))

    model_name = "Empirical" if args.model == 2 else "Empirical+F"
    print(f"[INFO] PAML AA model: {args.model} ({model_name})")
    print(f"[INFO] Amino-acid rate file: {args.aa_ratefile}")
    print(f"[INFO] Wrote ctl: {ctl_path}")
    print(f"[INFO] Running codeml in: {args.outdir}")

    try:
        subprocess.run(
            ["codeml", os.path.basename(ctl_path)],
            cwd=args.outdir,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        sys.exit(
            f"[ERROR] codeml failed with return code {e.returncode}"
        )

    print("[INFO] codeml finished successfully.")
    print(f"[INFO] Main output: {out_path}")
    print(f"[INFO] Ancestral reconstruction: {os.path.join(args.outdir, 'rst')}")


if __name__ == "__main__":
    main()
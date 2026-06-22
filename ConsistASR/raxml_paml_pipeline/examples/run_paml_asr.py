#!/usr/bin/env python3
"""
run_paml_asr.py

Minimal helper to:
  - generate a codeml .ctl file for amino-acid ASR
  - run codeml
  - keep all PAML outputs in a specified directory

This is intended as an *example* script for reproducing the
toy 7TM rhodopsin PAML ASR used in the ConsistASR examples.
"""

import argparse
import os
import subprocess
import textwrap
import sys


def main():
    ap = argparse.ArgumentParser(
        description="Generate a PAML codeml ctl file and run codeml for AA ASR."
    )
    ap.add_argument(
        "--aln",
        required=True,
        help="Alignment file in PHYLIP format (seqfile in codeml ctl).",
    )
    ap.add_argument(
        "--tree",
        required=True,
        help="Tree file in Newick format (treefile in codeml ctl).",
    )
    ap.add_argument(
        "--outdir",
        default="paml_asr",
        help="Output directory for ctl / codeml outputs (default: paml_asr).",
    )
    ap.add_argument(
        "--prefix",
        default="PAML_LG_FG",
        help="Prefix used for ctl / results files (default: PAML_LG_FG).",
    )
    ap.add_argument(
        "--aa-ratefile",
        default="lg.dat",
        help="Amino-acid rate file for model=2 (default: lg.dat).",
    )
    args = ap.parse_args()

    # basic checks
    for path in [args.aln, args.tree]:
        if not os.path.isfile(path):
            sys.exit(f"[ERROR] File not found: {path}")

    if not shutil_which("codeml"):
        sys.exit("[ERROR] codeml not found in PATH.")

    os.makedirs(args.outdir, exist_ok=True)
    ctl_path = os.path.join(args.outdir, f"{args.prefix}.ctl")
    out_path = os.path.join(args.outdir, f"{args.prefix}.out")

    # codeml works in the directory where the ctl is, with relative paths
    aln_rel = os.path.relpath(os.path.abspath(args.aln), args.outdir)
    tree_rel = os.path.relpath(os.path.abspath(args.tree), args.outdir)
    aa_rate_rel = os.path.relpath(os.path.abspath(args.aa_ratefile), args.outdir)

    ctl_text = f"""
    seqfile = {aln_rel}
    treefile = {tree_rel}
    outfile = {os.path.basename(out_path)}

    noisy = 9
    verbose = 1
    runmode = 0

    seqtype = 2     * 2: amino acid
    CodonFreq = 2

    model = 2       * Empirical (e.g. LG)
    aaRatefile = {aa_rate_rel}

    NSsites = 0
    icode = 0
    fix_alpha = 0
    alpha = 0.5
    ncatG = 4

    RateAncestor = 1   * output ancestral sequences
    cleandata = 1
    """

    with open(ctl_path, "w") as f:
        f.write(textwrap.dedent(ctl_text))

    print(f"[INFO] Wrote ctl: {ctl_path}")
    print(f"[INFO] Running codeml in: {args.outdir}")

    # run codeml with cwd=outdir
    try:
        subprocess.run(["codeml", os.path.basename(ctl_path)],
                       cwd=args.outdir,
                       check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"[ERROR] codeml failed with return code {e.returncode}")

    print("[INFO] codeml finished.")
    print(f"[INFO] Main output: {out_path}")


def shutil_which(cmd):
    """Small local 'which' wrapper (Python 3.6 compatibility)."""
    from shutil import which
    return which(cmd)


if __name__ == "__main__":
    main()

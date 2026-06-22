#!/usr/bin/env python3
"""
Create an ESR proxy alignment and tree for IQ-TREE ASR.

This script:
- reads an input gapped FASTA alignment and a Newick tree,
- removes one target extant sequence from the alignment,
- adds two dummy descendant sequences with the same gap pattern,
- replaces the target tip in the tree with a proxy subtree,
- labels the internal proxy node with the original target ID.

After running IQ-TREE ASR on the proxy dataset, the reconstructed state
for the internal node named TARGET can be compared with the original
held-out target sequence.
"""

import argparse
import re
from pathlib import Path


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)

        if name is not None:
            seqs[name] = "".join(buf)

    return seqs


def write_fasta(seqs: dict[str, str], path: Path, wrap: int = 80) -> None:
    with path.open("w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i + wrap] + "\n")


def validate_dummy_aa(value: str, label: str) -> str:
    aa = value.upper()
    if len(aa) != 1 or aa not in STANDARD_AA:
        raise SystemExit(
            f"[ERROR] {label} must be one standard amino-acid letter "
            f"({''.join(sorted(STANDARD_AA))})."
        )
    return aa


def make_dummy_sequence(target_seq: str, dummy_aa: str) -> str:
    return "".join("-" if c == "-" else dummy_aa for c in target_seq)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create an ESR proxy FASTA alignment and Newick tree by replacing "
            "one target tip with two dummy descendant tips."
        )
    )
    parser.add_argument("--msa", required=True,
                        help="Input gapped MSA in FASTA format.")
    parser.add_argument("--tree", required=True,
                        help="Input Newick tree with branch lengths.")
    parser.add_argument("--target", required=True,
                        help="Target tip ID to reconstruct; must match the FASTA header and tree tip label.")
    parser.add_argument("--dummy-bl", type=float, default=10.0,
                        help="Branch length for dummy leaves. Default: 10.0")
    parser.add_argument("--dummy-aa-a", default="L",
                        help="Dummy residue used for TARGET_A. Default: L")
    parser.add_argument("--dummy-aa-b", default="V",
                        help="Dummy residue used for TARGET_B. Default: V")
    parser.add_argument("--out-prefix", required=True,
                        help="Output prefix. Files written: <prefix>.fasta and <prefix>.tree")
    parser.add_argument("--wrap", type=int, default=80,
                        help="FASTA line width. Default: 80")
    args = parser.parse_args()

    msa_path = Path(args.msa)
    tree_path = Path(args.tree)
    target = args.target
    dummy_bl = args.dummy_bl
    dummy_aa_a = validate_dummy_aa(args.dummy_aa_a, "--dummy-aa-a")
    dummy_aa_b = validate_dummy_aa(args.dummy_aa_b, "--dummy-aa-b")

    if dummy_aa_a == dummy_aa_b:
        raise SystemExit("[ERROR] --dummy-aa-a and --dummy-aa-b should be different.")

    if args.wrap <= 0:
        raise SystemExit("[ERROR] --wrap must be a positive integer.")

    seqs = read_fasta(msa_path)
    if target not in seqs:
        raise SystemExit(f"[ERROR] Target '{target}' not found in MSA headers.")

    aln_len = len(next(iter(seqs.values())))
    if any(len(seq) != aln_len for seq in seqs.values()):
        raise SystemExit("[ERROR] MSA contains sequences with inconsistent lengths.")

    target_seq = seqs[target]
    dummy_seq_a = make_dummy_sequence(target_seq, dummy_aa_a)
    dummy_seq_b = make_dummy_sequence(target_seq, dummy_aa_b)

    a_name = f"{target}_A"
    b_name = f"{target}_B"

    new_seqs = {name: seq for name, seq in seqs.items() if name != target}
    if a_name in new_seqs or b_name in new_seqs:
        raise SystemExit("[ERROR] Dummy names already exist in MSA. Rename the target or existing sequences.")

    new_seqs[a_name] = dummy_seq_a
    new_seqs[b_name] = dummy_seq_b

    out_msa = Path(f"{args.out_prefix}.fasta")
    out_tree = Path(f"{args.out_prefix}.tree")

    write_fasta(new_seqs, out_msa, wrap=args.wrap)

    tree_str = tree_path.read_text().strip()

    if a_name in tree_str or b_name in tree_str:
        raise SystemExit("[ERROR] Dummy names already exist in the tree.")

    pattern = re.compile(
        rf"(?P<delim>[,(]){re.escape(target)}:(?P<brlen>[0-9eE.+-]+)"
    )
    matches = list(pattern.finditer(tree_str))

    if len(matches) == 0:
        raise SystemExit(
            f"[ERROR] Could not find '{target}:<branch_length>' in tree. "
            "Check that tree tip labels exactly match MSA headers and that the tree has branch lengths."
        )

    if len(matches) > 1:
        raise SystemExit(
            f"[ERROR] Found multiple occurrences of target '{target}' in tree. "
            "Tree tip labels should be unique."
        )

    def replace_target(match):
        return (
            f"{match.group('delim')}"
            f"({a_name}:{dummy_bl},{b_name}:{dummy_bl})"
            f"{target}:{match.group('brlen')}"
        )

    tree_new = pattern.sub(replace_target, tree_str, count=1)
    out_tree.write_text(tree_new + "\n")

    print("[DONE] Wrote:")
    print(f"  - {out_msa}")
    print(f"  - {out_tree}")
    print(
        f"[INFO] ESR internal node label = '{target}'; "
        f"dummy leaves = '{a_name}', '{b_name}'; "
        f"dummy BL = {dummy_bl}; dummy residues = {dummy_aa_a}/{dummy_aa_b}"
    )


if __name__ == "__main__":
    main()
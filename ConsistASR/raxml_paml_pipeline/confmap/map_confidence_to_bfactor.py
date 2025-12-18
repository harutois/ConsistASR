#!/usr/bin/env python3
import argparse
import sys
import math
import gemmi

# 3-letter to 1-letter amino acid code
AA3_TO_1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
    'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
    'TYR':'Y','VAL':'V','SEC':'U','PYL':'O'
}

def read_cif_sequence_and_plddt(cif_path, chain_id=None):
    """
    Read an AlphaFold mmCIF file and return:
      - structure sequence in 1-letter codes (for one polymer chain),
      - a list of lists of atoms per residue (for that chain),
      - the gemmi.Structure object,
      - the actual chain ID used.

    If chain_id is None, the longest polymer chain is selected.
    """
    st = gemmi.read_structure(cif_path)
    mdl = st[0]
    chains = [ch for ch in mdl]

    def is_polymer_chain(ch):
        for res in ch:
            if getattr(res, "entity_type", None) == gemmi.EntityType.Polymer:
                return True
        return False

    # --- chain selection ---
    if chain_id:
        ch = next((c for c in chains if c.name.strip() == chain_id.strip()), None)
        if ch is None:
            sys.exit(f"[ERROR] Chain '{chain_id}' not found in structure.")
        if not is_polymer_chain(ch):
            sys.exit(f"[ERROR] Chain '{chain_id}' is not a polymer chain.")
    else:
        polys = [c for c in chains if is_polymer_chain(c)]
        if not polys:
            sys.exit("[ERROR] No polymer chain found in structure.")
        # pick the polymer chain with the largest number of residues
        ch = max(polys, key=lambda c: sum(1 for _ in c))
        chain_id = ch.name

    # --- build sequence and atom list per residue ---
    seq = []
    res_atoms = []
    for res in ch:
        if getattr(res, "entity_type", None) != gemmi.EntityType.Polymer:
            continue
        aa3 = res.name.strip().upper()
        aa1 = AA3_TO_1.get(aa3, 'X')
        seq.append(aa1)
        atoms = [atom for atom in res]
        res_atoms.append(atoms)

    return ''.join(seq), res_atoms, st, chain_id

def read_fasta_first_seq(path):
    """Read the first sequence (concatenated) from a FASTA file."""
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return ''.join(seq)

def read_pp_values(path):
    """Read per-site PP values (one float per line)."""
    vals = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            try:
                vals.append(float(s))
            except ValueError:
                sys.exit(f"[ERROR] Non-numeric PP value in {path}: {s}")
    return vals

def make_ungapped_to_aln_map(aln_seq):
    """
    From an aligned sequence with gaps ('-'), return:
      - ungapped sequence (string),
      - a list mapping ungapped indices (0-based) to
        alignment columns (0-based).
    """
    map_list = []
    ungapped = []
    for i, c in enumerate(aln_seq):
        if c != '-':
            map_list.append(i)
            ungapped.append(c)
    return ''.join(ungapped), map_list

def report_seq_mismatch(name_a, seq_a, name_b, seq_b, max_show=50):
    """Print a brief report if two sequences do not match."""
    if seq_a == seq_b:
        return
    print(f"[WARN] Sequence mismatch between {name_a} and {name_b}.")
    diffs = []
    for i, (x, y) in enumerate(zip(seq_a, seq_b)):
        if x != y:
            diffs.append((i, x, y))
            if len(diffs) >= max_show:
                break
    msg = ', '.join([f"{i}:{x}->{y}" for i, x, y in diffs]) if diffs else "length or trailing mismatch"
    print(f"       First diffs: {msg}")
    if len(seq_a) != len(seq_b):
        print(f"       Lengths: {name_a}={len(seq_a)}, {name_b}={len(seq_b)}")

def scale_pp_to_100(pp_vals):
    """
    Automatically scale PP values to 0–100:
      - if max(pp) <= 1.0, treat as 0–1 and multiply by 100,
      - otherwise assume values are already in 0–100.
    """
    mx = max(pp_vals) if pp_vals else 1.0
    if mx <= 1.00001:   # assume 0–1
        return [v*100.0 for v in pp_vals], "0–1→×100"
    else:
        return pp_vals, "already 0–100"

def write_pdb_with_b(st, res_atoms, bvals, out_path):
    """
    Write a PDB file where all atoms in residue i get
    the same B-factor bvals[i].
    """
    idx = 0
    for atoms in res_atoms:
        if idx >= len(bvals):
            print(f"[WARN] B-factor array shorter than residues; skipping from index {idx}")
            break
        bval = bvals[idx]
        for atom in atoms:
            atom.b_iso = bval
        idx += 1

    with open(out_path, "w") as f:
        f.write(st.make_pdb_string())
    print(f"[OK] Wrote: {out_path}")

def main():
    ap = argparse.ArgumentParser(
        description="Map ASR site-wise PP and/or AlphaFold pLDDT to B-factor, handling alignment gaps."
    )
    ap.add_argument('--cif', required=True,
                    help='AlphaFold mmCIF file')
    ap.add_argument('--asr_fasta', required=True,
                    help='ASR sequence in FASTA (aligned, with gaps \"-\"; first sequence is used)')
    ap.add_argument('--pp_file', required=True,
                    help='Per-site PP values file (one value per alignment column)')
    ap.add_argument('--chain', default=None,
                    help='Chain ID in structure (default: longest polymer chain)')
    ap.add_argument('--mode',
                    choices=['plddt', 'pp', 'diff', 'prod'],
                    required=True,
                    help=('Which per-residue metric to map to B-factor: '
                          'plddt (=original AlphaFold pLDDT), '
                          'pp (=ASR posterior probability), '
                          'diff (=pLDDT−PP), '
                          'prod (=PP×pLDDT/100)'))
    ap.add_argument('--out', required=True,
                    help='Output PDB path')
    args = ap.parse_args()

    # 1) Read structure and sequence
    struct_seq, res_atoms, st, chain_id = read_cif_sequence_and_plddt(args.cif, args.chain)
    print(f"[INFO] Structure chain={chain_id}, residues={len(res_atoms)}")

    # 2) Read ASR alignment and PP
    asr_aln = read_fasta_first_seq(args.asr_fasta).strip()
    if not asr_aln:
        sys.exit("[ERROR] Empty ASR fasta.")
    pp_vals = read_pp_values(args.pp_file)
    if len(pp_vals) != len(asr_aln):
        print(f"[WARN] PP length ({len(pp_vals)}) != alignment length ({len(asr_aln)}). Truncating to min.")
        L = min(len(pp_vals), len(asr_aln))
        asr_aln = asr_aln[:L]
        pp_vals = pp_vals[:L]

    # 3) Build mapping from ungapped ASR sequence to alignment columns
    asr_ungap, map_u2aln = make_ungapped_to_aln_map(asr_aln)
    print(f"[INFO] Alignment length={len(asr_aln)}, ungapped length={len(asr_ungap)}")

    # 4) Compare ASR(ungapped) vs structure sequence
    report_seq_mismatch("ASR(ungapped)", asr_ungap, "Structure", struct_seq)
    if len(asr_ungap) != len(struct_seq):
        print("[ERROR] Ungapped ASR length != structure sequence length; cannot safely map B-factors.")
        sys.exit(1)

    # 5) Re-map PP (0–100) to per-residue order
    pp_scaled, scale_info = scale_pp_to_100(pp_vals)
    print(f"[INFO] PP scale: {scale_info}")
    per_res_pp = []
    for u, aln_col in enumerate(map_u2aln):
        v = pp_scaled[aln_col]
        v = 0.0 if math.isnan(v) else max(0.0, min(100.0, v))
        per_res_pp.append(v)

    # 6) Compute per-residue pLDDT (mean B-factor per residue)
    per_res_plddt = []
    for atoms in res_atoms:
        if not atoms:
            per_res_plddt.append(0.0)
            continue
        mean_b = sum(a.b_iso for a in atoms) / len(atoms)
        per_res_plddt.append(mean_b)

    # 7) Decide which metric to map to B-factor
    if args.mode == 'pp':
        bvals = per_res_pp
    elif args.mode == 'plddt':
        bvals = per_res_plddt
    elif args.mode == 'diff':
        # Here we store pLDDT - PP (both assumed 0–100)
        # Any scaling to 0–100 can be applied externally if desired.
        bvals = [pl - pp for pl, pp in zip(per_res_plddt, per_res_pp)]
    else:  # prod
        # PP × pLDDT / 100 → 0–100
        bvals = [(pp * pl) / 100.0 for pp, pl in zip(per_res_pp, per_res_plddt)]

    # 8) Write PDB
    write_pdb_with_b(st, res_atoms, bvals, args.out)

if __name__ == '__main__':
    main()

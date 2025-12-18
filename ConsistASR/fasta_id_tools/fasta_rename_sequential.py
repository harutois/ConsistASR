#!/usr/bin/env python3
import argparse
import sys
from textwrap import wrap

def read_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        yield header, ''.join(seq_chunks)

def write_fasta(records, out_path):
    with open(out_path, 'w') as out:
        for h, s in records:
            out.write(f">{h}\n")
            for chunk in wrap(s, 60):
                out.write(chunk + "\n")

def main():
    ap = argparse.ArgumentParser(
        description="Rename all FASTA IDs to seq01, seq02, ... and output a mapping table."
    )
    ap.add_argument("--in", dest="in_fasta", required=True, help="Input FASTA")
    ap.add_argument("--out", dest="out_fasta", required=True, help="Output FASTA with new IDs")
    ap.add_argument("--map", dest="map_tsv", required=True, help="Output mapping TSV (new_id<TAB>old_header)")
    args = ap.parse_args()

    records = list(read_fasta(args.in_fasta))
    n = len(records)
    if n == 0:
        sys.exit("[ERROR] No sequences found in input FASTA.")

    # Digit Determination (Example: 1–99 -> 2 digits, 100–999 -> 3 digits)
    width = max(2, len(str(n)))

    new_records = []
    mapping_lines = []
    used_ids = set()

    for i, (old_header, seq) in enumerate(records, start=1):
        base_id = f"seq{str(i).zfill(width)}"
        new_id = base_id
        suffix = 1
        while new_id in used_ids:
            new_id = f"{base_id}_{suffix}"
            suffix += 1
        used_ids.add(new_id)

        new_records.append((new_id, seq))
        mapping_lines.append(f"{new_id}\t{old_header}")

    write_fasta(new_records, args.out_fasta)
    with open(args.map_tsv, "w") as m:
        m.write("new_id\told_header\n")
        for line in mapping_lines:
            m.write(line + "\n")

    print(f"[INFO] Wrote {len(new_records)} sequences to {args.out_fasta}")
    print(f"[INFO] Wrote mapping table to {args.map_tsv}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import sys
import re
from textwrap import wrap

VALID_RE = re.compile(r'^([A-Za-z0-9_]+)')

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
        description="Truncate FASTA IDs at the first non-alnum/underscore character and output mapping."
    )
    ap.add_argument("--in", dest="in_fasta", required=True, help="Input FASTA")
    ap.add_argument("--out", dest="out_fasta", required=True, help="Output FASTA with truncated IDs")
    ap.add_argument("--map", dest="map_tsv", required=True, help="Output mapping TSV (new_id<TAB>old_header)")
    args = ap.parse_args()

    records = list(read_fasta(args.in_fasta))
    if not records:
        sys.exit("[ERROR] No sequences found in input FASTA.")

    new_records = []
    mapping_lines = []
    used_ids = set()
    auto_id_counter = 1

    for old_header, seq in records:
        m = VALID_RE.match(old_header)
        if m:
            base = m.group(1)
        else:
            # Fallback for abnormal cases where the start does not begin with alphanumeric characters or an underscore
            base = f"id{auto_id_counter}"
            auto_id_counter += 1
            print(f"[WARN] Header '{old_header}' has no leading [A-Za-z0-9_]; using fallback ID '{base}'", file=sys.stderr)

        new_id = base
        suffix = 1
        while new_id in used_ids:
            # When multiple instances of the same base appear, use _1, _2,... to make them unique
            new_id = f"{base}_{suffix}"
            suffix += 1
        if new_id != base:
            print(f"[WARN] Duplicate truncated ID '{base}' detected, renamed to '{new_id}'", file=sys.stderr)

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

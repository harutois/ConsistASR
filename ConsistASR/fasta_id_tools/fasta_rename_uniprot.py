#!/usr/bin/env python3
import argparse
import sys
import re
from textwrap import wrap

UNIPROT_RE = re.compile(
    r'^(?P<db>sp|tr)\|(?P<accession>[^|]+)\|(?P<entry_name>\S+)\s*(?P<rest>.*)$'
)

FIELD_MARKER_RE = re.compile(r'\s(?:OS|OX|GN|PE|SV)=')

SAFE_ID_RE = re.compile(r'[^A-Za-z0-9_]+')


def read_fasta(path):
    header = None
    seq_chunks = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                if header is None:
                    sys.exit("[ERROR] Sequence line found before the first FASTA header.")
                seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)


def write_fasta(records, out_path, width=60):
    with open(out_path, "w", encoding="utf-8") as out:
        for h, s in records:
            out.write(f">{h}\n")
            for chunk in wrap(s, width):
                out.write(chunk + "\n")


def clean_id(text):
    text = SAFE_ID_RE.sub("_", text)
    text = text.strip("_")
    return text if text else "unknown"


def clean_tsv_field(text):
    if text is None:
        return ""
    return str(text).replace("\t", " ").replace("\n", " ").strip()


def positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("Value must be a positive integer.")
    return ivalue


def extract_uniprot_field(rest, tag):
    """
    Extract fields such as OS=..., OX=..., GN=..., PE=..., SV=...
    from the UniProt FASTA description.
    """
    pattern = re.compile(rf'(?:^|\s){tag}=([^=]*?)(?=\s(?:OS|OX|GN|PE|SV)=|$)')
    m = pattern.search(rest)
    if not m:
        return ""
    return m.group(1).strip()


def extract_protein_name(rest):
    """
    UniProt FASTA format:
    sp|ACCESSION|ENTRY_NAME Protein name OS=Organism OX=TaxID GN=Gene PE=... SV=...
    Protein name is the part before the first UniProt key-value marker.
    """
    m = FIELD_MARKER_RE.search(rest)
    if not m:
        return rest.strip()
    return rest[:m.start()].strip()


def parse_uniprot_header(header):
    m = UNIPROT_RE.match(header)
    if not m:
        return None

    db = m.group("db")
    accession = m.group("accession")
    entry_name = m.group("entry_name")
    rest = m.group("rest").strip()

    protein_name = extract_protein_name(rest)
    organism = extract_uniprot_field(rest, "OS")
    ox = extract_uniprot_field(rest, "OX")
    gene = extract_uniprot_field(rest, "GN")
    pe = extract_uniprot_field(rest, "PE")
    sv = extract_uniprot_field(rest, "SV")

    return {
        "db": db,
        "accession": accession,
        "entry_name": entry_name,
        "protein_name": protein_name,
        "organism": organism,
        "ox": ox,
        "gene": gene,
        "pe": pe,
        "sv": sv,
    }


def make_base_id(info, id_format, separator):
    if id_format == "accession":
        base = info["accession"]
    elif id_format in ("entry", "entry_name"):
        base = info["entry_name"]
    elif id_format == "accession_entry":
        base = f'{info["accession"]}{separator}{info["entry_name"]}'
    elif id_format == "db_accession_entry":
        base = f'{info["db"]}{separator}{info["accession"]}{separator}{info["entry_name"]}'
    else:
        raise ValueError(f"Unsupported id format: {id_format}")

    return clean_id(base)


def make_fallback_id(header, fallback_mode, counter):
    if fallback_mode == "sanitize":
        return clean_id(header), counter
    if fallback_mode == "sequential":
        fallback_id = f"id{counter:05d}"
        return fallback_id, counter + 1
    if fallback_mode == "strict":
        sys.exit(f"[ERROR] Non-UniProt FASTA header found: {header}")
    raise ValueError(f"Unsupported fallback mode: {fallback_mode}")


def uniquify_id(base, used_ids):
    new_id = base
    suffix = 1
    duplicate_resolved = "no"

    while new_id in used_ids:
        new_id = f"{base}_{suffix}"
        suffix += 1
        duplicate_resolved = "yes"

    used_ids.add(new_id)
    return new_id, duplicate_resolved


def main():
    ap = argparse.ArgumentParser(
        description=(
    "Rename UniProt FASTA headers to short IDs and output a mapping table. "
    "Default ID format: ENTRY_NAME."
    )
    )
    ap.add_argument("--in", dest="in_fasta", required=True, help="Input FASTA")
    ap.add_argument("--out", dest="out_fasta", required=True, help="Output FASTA with renamed IDs")
    ap.add_argument("--map", dest="map_tsv", required=True, help="Output mapping TSV")
    ap.add_argument(
    "--id-format",
    choices=["entry", "entry_name", "accession", "accession_entry", "db_accession_entry"],
    default="entry",
    help=(
        "New FASTA ID format. "
        "Default: entry. "
        "Use accession_entry to include UniProt accession."
    ),
    )
    ap.add_argument(
        "--separator",
        default="_",
        help="Separator used for composite IDs. Default: '_'",
    )
    ap.add_argument(
        "--prefix",
        default="",
        help="Optional prefix added to every new ID, e.g. qNOR_",
    )
    ap.add_argument(
        "--fallback",
        choices=["sanitize", "sequential", "strict"],
        default="sanitize",
        help=(
            "How to handle non-UniProt headers. "
            "sanitize: replace unsafe characters; sequential: id00001; strict: stop with error. "
            "Default: sanitize"
        ),
    )
    ap.add_argument(
    "--wrap",
    type=positive_int,
    default=60,
    help="Line width for output FASTA sequences. Default: 60",
    )
    args = ap.parse_args()

    records = list(read_fasta(args.in_fasta))
    if not records:
        sys.exit("[ERROR] No sequences found in input FASTA.")

    new_records = []
    mapping_rows = []
    used_ids = set()
    fallback_counter = 1

    n_uniprot = 0
    n_fallback = 0
    n_duplicate = 0

    for old_header, seq in records:
        info = parse_uniprot_header(old_header)

        if info is not None:
            n_uniprot += 1
            base = make_base_id(info, args.id_format, args.separator)
            fallback_used = "no"
        else:
            n_fallback += 1
            base, fallback_counter = make_fallback_id(old_header, args.fallback, fallback_counter)
            info = {
                "db": "",
                "accession": "",
                "entry_name": "",
                "protein_name": "",
                "organism": "",
                "ox": "",
                "gene": "",
                "pe": "",
                "sv": "",
            }
            fallback_used = args.fallback
            print(
                f"[WARN] Non-UniProt header handled by fallback='{args.fallback}': {old_header}",
                file=sys.stderr,
            )

        if args.prefix:
            base = clean_id(f"{args.prefix}{base}")

        new_id, duplicate_resolved = uniquify_id(base, used_ids)
        if duplicate_resolved == "yes":
            n_duplicate += 1
            print(
                f"[WARN] Duplicate ID base '{base}' detected; renamed to '{new_id}'",
                file=sys.stderr,
            )

        new_records.append((new_id, seq))

        mapping_rows.append([
            new_id,
            old_header,
            info["db"],
            info["accession"],
            info["entry_name"],
            info["protein_name"],
            info["organism"],
            info["ox"],
            info["gene"],
            info["pe"],
            info["sv"],
            fallback_used,
            duplicate_resolved,
        ])

    write_fasta(new_records, args.out_fasta, width=args.wrap)

    with open(args.map_tsv, "w", encoding="utf-8") as m:
        m.write(
            "\t".join([
                "new_id",
                "old_header",
                "db",
                "accession",
                "entry_name",
                "protein_name",
                "organism",
                "ox",
                "gene",
                "pe",
                "sv",
                "fallback_used",
                "duplicate_resolved",
            ]) + "\n"
        )
        for row in mapping_rows:
            m.write("\t".join(clean_tsv_field(x) for x in row) + "\n")

    print(f"[INFO] Input sequences: {len(records)}")
    print(f"[INFO] UniProt-format headers: {n_uniprot}")
    print(f"[INFO] Fallback-renamed headers: {n_fallback}")
    print(f"[INFO] Duplicate IDs resolved: {n_duplicate}")
    print(f"[INFO] Wrote renamed FASTA to: {args.out_fasta}")
    print(f"[INFO] Wrote mapping table to: {args.map_tsv}")


if __name__ == "__main__":
    main()
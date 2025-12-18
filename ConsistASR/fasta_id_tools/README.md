# FASTA ID utilities

This directory provides small helper scripts to **normalize FASTA sequence IDs** before passing alignments to tools that are sensitive to spaces or special characters in sequence names (e.g. PHYLIP, PAML, RAxML).

All scripts:

- take a multi-FASTA file as input
- write a new FASTA with modified IDs
- write a tab-separated mapping file (`new_id<TAB>old_header`) so that original IDs can be recovered if needed
- ensure that new IDs are **unique** (adding `_1`, `_2`, … if necessary)

Typical use case:  
cleaning FASTA headers such as

```text
>AM_5S_00009|AM_5S-c1__COG5524_Bacteriorhodopsin_bac_rhodopsin_ASGARD
>QBQ84358.1 schizorhodopsin [Promethearchaeota archaeon]
>SAMEA 2622822_312577
```

before converting to PHYLIP / feeding into ASR pipelines.

---

## 1. `fasta_rename_sequential.py`

Rename **all** sequences to `seq01`, `seq02`, … and write a mapping table.

This is the safest and most PHYLIP-friendly option when you do not need human-readable IDs inside the ASR tools.

### Usage

```bash
python3 fasta_rename_sequential.py \
  --in  input.fasta \
  --out renamed.fasta \
  --map id_map.tsv
```

* `--in` : input multi-FASTA
* `--out`: output FASTA with IDs like `seq01`, `seq02`, …
* `--map`: mapping table (TSV) with header line

Example of `id_map.tsv`:

```text
new_id  old_header
seq01   AM_5S_00009|AM_5S-c1__COG5524_Bacteriorhodopsin_bac_rhodopsin_ASGARD
seq02   QBQ84358.1 schizorhodopsin [Promethearchaeota archaeon]
seq03   SAMEA 2622822_312577
```

The number of digits in `seqXX` is automatically chosen based on the number of sequences (e.g. `seq01`–`seq99`, `seq001`–…).

If a generated ID accidentally collides (should be extremely rare), `_1`, `_2`, … are appended and a warning is printed to STDERR.

---

## 2. `fasta_truncate_id_simple.py`

Truncate each header at the **first character that is not** `[A–Z, a–z, 0–9, _]`, and use the remaining prefix as the new ID.

This keeps IDs relatively short while preserving a recognizable fragment of the original identifier.
Examples:

* `AM_5S_00009|AM_5S-c1__COG5524_Bacteriorhodopsin_bac_rhodopsin_ASGARD` → `AM_5S_00009`
* `QBQ84358.1 schizorhodopsin [Promethearchaeota archaeon]` → `QBQ84358`
* `SAMEA 2622822_312577` → `SAMEA`

### Usage

```bash
python3 fasta_truncate_id_simple.py \
  --in  input.fasta \
  --out truncated.fasta \
  --map id_map_truncate.tsv
```

* `--in` : input multi-FASTA
* `--out`: FASTA with truncated IDs (single token, no spaces)
* `--map`: mapping table (`new_id<TAB>old_header`)

If two different headers would lead to the same truncated ID (e.g. `QBQ84358.1` and `QBQ84358.2` → both `QBQ84358`), the script automatically generates `QBQ84358_1`, `QBQ84358_2`, … and prints a warning to STDERR.

If a header does **not** start with `[A–Z0–9_]` at all (pathological case), a fallback ID (`id1`, `id2`, …) is used and a warning is emitted.

---

## 3. `fasta_sanitize_id_underscore.py`

Replace every non-`[A–Z, a–z, 0–9, _]` character by an underscore, compress multiple underscores, and strip leading/trailing underscores.

This preserves most of the original header information, but makes it safe for programs that dislike spaces or punctuation.

Examples:

* `AM_5S_00009|AM_5S-c1__COG5524_Bacteriorhodopsin_bac_rhodopsin_ASGARD`
  → `AM_5S_00009_AM_5S_c1__COG5524_Bacteriorhodopsin_bac_rhodopsin_ASGARD`
* `QBQ84358.1 schizorhodopsin [Promethearchaeota archaeon]`
  → `QBQ84358_1_schizorhodopsin__Promethearchaeota_archaeon_`
* `SAMEA 2622822_312577` → `SAMEA_2622822_312577`

Note that resulting IDs can be quite long; for strict PHYLIP use you may prefer script (1) or (2).

### Usage

```bash
python3 fasta_sanitize_id_underscore.py \
  --in  input.fasta \
  --out sanitized.fasta \
  --map id_map_sanitize.tsv
```

If the sanitized ID becomes empty (e.g. header was only symbols), the script falls back to `id1`, `id2`, … and warns on STDERR.
If duplicates are produced, `_1`, `_2`, … are appended with warnings.

---

## 4. Recommended usage

Typical workflow options:

* **Maximum safety for PHYLIP / PAML**
  Use **sequential renaming**:

  ```bash
  python3 fasta_rename_sequential.py \
    --in  raw.fasta \
    --out phylip_safe.fasta \
    --map phylip_safe_map.tsv
  ```

* **Keep short but meaningful IDs** (good compromise for this ConsistASR project)
  Use **truncation**:

  ```bash
  python3 fasta_truncate_id_simple.py \
    --in  raw.fasta \
    --out trimmed_ids.fasta \
    --map trimmed_ids_map.tsv
  ```

* **Debugging / manual inspection** with “human-readable” IDs that still avoid spaces and punctuation
  Use **underscore sanitization**:

  ```bash
  python3 fasta_sanitize_id_underscore.py \
    --in  raw.fasta \
    --out underscored_ids.fasta \
    --map underscored_ids_map.tsv
  ```

You can also chain these scripts, for example:

1. First run `fasta_sanitize_id_underscore.py` to normalize headers.
2. Then run `fasta_rename_sequential.py` if you need very short IDs for PHYLIP, while keeping the two mapping files as documentation.

---

## 5. Notes

* All scripts assume standard multi-FASTA format (`>header` lines followed by one or more sequence lines).
* Sequences are wrapped at 60 characters per line in the output.
* Mapping files always include a header line: `new_id<TAB>old_header`.
* Warnings related to duplicated IDs or unusual headers are written to STDERR, so they will appear in the terminal but not in the mapping file itself.

---

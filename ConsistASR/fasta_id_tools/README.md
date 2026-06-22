# FASTA ID utilities

This directory provides small helper scripts to **normalize FASTA sequence IDs** before passing alignments to tools that are sensitive to spaces or special characters in sequence names (e.g. PHYLIP, PAML, RAxML).

All scripts:

- take a multi-FASTA file as input
- write a new FASTA with modified IDs
- write a tab-separated mapping file so that original IDs can be recovered if needed
  - the simple scripts write `new_id<TAB>old_header`
  - `fasta_rename_uniprot.py` writes an extended UniProt-aware mapping table
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

## 4. `fasta_rename_uniprot.py`

Rename UniProt-style FASTA headers to short, tool-friendly IDs and write an extended mapping table.

This script is intended for FASTA files downloaded from UniProtKB, where headers typically look like:

```text
>tr|A0A285P2K2|A0A285P2K2_NATPI Nitric oxide reductase, NorZ apoprotein OS=Natronoarchaeum philippinense OX=558529 GN=SAMN06269185_2558 PE=4 SV=1
>tr|A0A871BJ59|A0A871BJ59_HALGI Homolog to nitric oxide reductase subunit B OS=Haloferax gibbonsii OX=35746 GN=HfgLR_13555 PE=4 SV=1
>tr|A0ABD5PKQ0|A0ABD5PKQ0_9EURY Nitric-oxide reductase large subunit OS=Halosolutus amylolyticus OX=2932267 GN=ACFO5R_03125 PE=4 SV=1
>tr|B9LQ42|B9LQ42_HALLT Cytochrome b subunit of nitric oxide reductase OS=Halorubrum lacusprofundi (strain ATCC 49239 / DSM 5036 / JCM 8891 / ACAM 34) OX=416348 GN=Hlac_1902 PE=4 SV=1
```

By default, the script uses the UniProt **entry name** as the new FASTA ID:

```text
tr|A0A285P2K2|A0A285P2K2_NATPI ... → A0A285P2K2_NATPI
tr|A0A871BJ59|A0A871BJ59_HALGI ... → A0A871BJ59_HALGI
tr|A0ABD5PKQ0|A0ABD5PKQ0_9EURY ... → A0ABD5PKQ0_9EURY
tr|B9LQ42|B9LQ42_HALLT ...         → B9LQ42_HALLT
```

This is usually short enough for tree visualization and still informative for manual inspection.

### Usage

```bash
python3 fasta_rename_uniprot.py \
  --in  uniprot_input.fasta \
  --out uniprot_short.fasta \
  --map uniprot_id_map.tsv
```

Default output ID format:

```text
ENTRY_NAME
```

For example:

```text
A0A285P2K2_NATPI
A0A871BJ59_HALGI
A0ABD5PKQ0_9EURY
B9LQ42_HALLT
```

### Alternative ID formats

If you want to include the UniProt accession in the FASTA ID, use:

```bash
python3 fasta_rename_uniprot.py \
  --in  uniprot_input.fasta \
  --out uniprot_accession_entry.fasta \
  --map uniprot_accession_entry_map.tsv \
  --id-format accession_entry
```

This gives IDs such as:

```text
A0A285P2K2_A0A285P2K2_NATPI
A0A871BJ59_A0A871BJ59_HALGI
A0ABD5PKQ0_A0ABD5PKQ0_9EURY
B9LQ42_B9LQ42_HALLT
```

Available `--id-format` options:

```text
entry              default; use UniProt entry name only
entry_name         alias of entry
accession          use UniProt accession only
accession_entry    use ACCESSION_ENTRY_NAME
db_accession_entry use sp/tr_ACCESSION_ENTRY_NAME
```

### Optional prefix for UniProt-derived subsets

`fasta_rename_uniprot.py` also supports an optional `--prefix` argument.
This option is specific to the UniProt-aware renaming script and is useful when separately collected UniProt-derived subsets may later be combined.

For example, qNOR and cNOR candidate sequences can be given different prefixes:

```bash
python3 fasta_rename_uniprot.py \
  --in  uniprot_qNOR.fasta \
  --out qNOR_short.fasta \
  --map qNOR_id_map.tsv \
  --prefix qNOR_
```

Example:

```text
A0A285P2K2_NATPI → qNOR_A0A285P2K2_NATPI
```

The prefix is sanitized together with the generated ID, so the final IDs contain only letters, numbers, and underscores.


### Mapping table

Unlike the other simple FASTA ID utilities, this script writes an extended mapping table.

Example columns:

```text
new_id
old_header
db
accession
entry_name
protein_name
organism
ox
gene
pe
sv
fallback_used
duplicate_resolved
```

This makes it easier to recover the original UniProt accession, organism name, gene name, and protein annotation after sequence renaming.

Example:

```text
new_id	old_header	db	accession	entry_name	protein_name	organism	ox	gene	pe	sv	fallback_used	duplicate_resolved
A0A285P2K2_NATPI	tr|A0A285P2K2|A0A285P2K2_NATPI Nitric oxide reductase, NorZ apoprotein OS=Natronoarchaeum philippinense OX=558529 GN=SAMN06269185_2558 PE=4 SV=1	tr	A0A285P2K2	A0A285P2K2_NATPI	Nitric oxide reductase, NorZ apoprotein	Natronoarchaeum philippinense	558529	SAMN06269185_2558	4	1	no	no
```

### Handling duplicate IDs

If two headers generate the same new ID, the script appends `_1`, `_2`, … and prints a warning to STDERR.

Example:

```text
A0A285P2K2_NATPI
A0A285P2K2_NATPI_1
A0A285P2K2_NATPI_2
```

The mapping table records whether duplicate resolution was used in the `duplicate_resolved` column.

### Handling non-UniProt headers

By default, non-UniProt headers are sanitized rather than causing the script to stop.

Available fallback modes:

```text
--fallback sanitize     default; replace unsafe characters with underscores
--fallback sequential   use id00001, id00002, ...
--fallback strict       stop with an error if a non-UniProt header is found
```

For strict checking of UniProt-only FASTA files:

```bash
python3 fasta_rename_uniprot.py \
  --in  uniprot_input.fasta \
  --out uniprot_short.fasta \
  --map uniprot_id_map.tsv \
  --fallback strict
```

This is useful when you want to ensure that all input headers follow the expected UniProt format.

Small example input and expected output files are provided in:

```text
fasta_id_tools/examples/uniprot_rename/
```

---

## 5. Recommended usage

Typical workflow options:

* **Maximum safety for PHYLIP / PAML**
  Use **sequential renaming**:

  ```bash
  python3 fasta_rename_sequential.py \
    --in  raw.fasta \
    --out phylip_safe.fasta \
    --map phylip_safe_map.tsv
  ```

* **Keep short but meaningful IDs for general FASTA files**
  Use **simple truncation**:

  ```bash
  python3 fasta_truncate_id_simple.py \
    --in  raw.fasta \
    --out trimmed_ids.fasta \
    --map trimmed_ids_map.tsv
  ```

* **Keep UniProt entry names while preserving accession and annotation information**
  Use **UniProt-aware renaming**:

  ```bash
  python3 fasta_rename_uniprot.py \
    --in  uniprot_raw.fasta \
    --out uniprot_short.fasta \
    --map uniprot_id_map.tsv
  ```

  This is recommended for UniProtKB FASTA files used in ASR pipelines.
  When multiple UniProt-derived subsets are renamed separately and later merged, `--prefix` can be used to label each subset, for example `qNOR_` or `cNOR_`.

* **Debugging / manual inspection** with “human-readable” IDs that still avoid spaces and punctuation
  Use **underscore sanitization**:

  ```bash
  python3 fasta_sanitize_id_underscore.py \
    --in  raw.fasta \
    --out underscored_ids.fasta \
    --map underscored_ids_map.tsv
  ```

You can also chain these scripts, for example:

1. First run `fasta_rename_uniprot.py` to shorten UniProt headers while keeping an extended mapping table.
2. Then run `fasta_rename_sequential.py` if you need very short IDs for strict PHYLIP/PAML workflows.
3. Keep both mapping files as documentation.

---

## 6. Notes

* All scripts assume standard multi-FASTA format (`>header` lines followed by one or more sequence lines).

* Sequences are wrapped at 60 characters per line in the output.

* All scripts ensure that new IDs are unique by appending `_1`, `_2`, … if necessary.

* The simple scripts write mapping files with at least:

  ```text
  new_id<TAB>old_header
  ```

* `fasta_rename_uniprot.py` writes an extended mapping table containing UniProt-specific fields such as accession, entry name, organism, gene name, protein name, PE, and SV.

* Warnings related to duplicated IDs, unusual headers, or fallback renaming are written to STDERR, so they will appear in the terminal but not interrupt the mapping file output unless `--fallback strict` is used.

---

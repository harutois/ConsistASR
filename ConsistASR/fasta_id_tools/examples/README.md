# Example: UniProt FASTA renaming

This directory contains a small example dataset for `fasta_rename_uniprot.py`.

The input FASTA file contains UniProt-style headers:

```text
input_uniprot.fasta
```

Run the script from the `fasta_id_tools` directory:

```bash
python3 fasta_rename_uniprot.py \
  --in  examples/uniprot_rename/input_uniprot.fasta \
  --out examples/uniprot_rename/output_entry.fasta \
  --map examples/uniprot_rename/output_entry_map.tsv
```

By default, UniProt entry names are used as the new FASTA IDs.

The expected output files are:

```text
expected_entry.fasta
expected_entry_map.tsv
```

You can compare the generated files with the expected files using:

```bash
diff examples/uniprot_rename/output_entry.fasta examples/uniprot_rename/expected_entry.fasta
diff examples/uniprot_rename/output_entry_map.tsv examples/uniprot_rename/expected_entry_map.tsv
```

If the script behaves as expected, `diff` should produce no output.

The generated files can be removed after testing:

```bash
rm examples/uniprot_rename/output_entry.fasta
rm examples/uniprot_rename/output_entry_map.tsv
```

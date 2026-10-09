# Example: UniProt FASTA renaming

This directory contains a small example dataset and expected outputs for testing `fasta_rename_uniprot.py`.

The input FASTA file contains UniProt-style headers:

```text
input_uniprot.fasta
```

Run the script from the `fasta_id_tools` directory:

```bash
python3 fasta_rename_uniprot.py \
  --in  examples/input_uniprot.fasta \
  --out examples/output_entry.fasta \
  --map examples/output_entry_map.tsv
```

By default, UniProt entry names are used as the new FASTA IDs.

The corresponding reference outputs provided in this directory are:

```text
expected_entry.fasta
expected_entry_map.tsv
```

You can compare the generated files with the reference outputs using:

```bash
diff examples/output_entry.fasta examples/expected_entry.fasta
diff examples/output_entry_map.tsv examples/expected_entry_map.tsv
```

If the script behaves as expected, `diff` should produce no output.

The generated files can be removed after testing:

```bash
rm examples/output_entry.fasta
rm examples/output_entry_map.tsv
```

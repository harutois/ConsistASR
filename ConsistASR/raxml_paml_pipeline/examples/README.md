## Toy example: 7TM rhodopsin dataset

This directory contains a small 7TM rhodopsin toy dataset and bundled example outputs for the RAxML + PAML-based ConsistASR workflow.

* `toy_7tm_rhodopsin.fasta`

Example outputs (trees, ASR state files, indel-aware FASTA, confmap PDBs) are provided under the corresponding `examples/` subdirectories, so users can:

1. Inspect intermediate files
2. Re-run the scripts step by step
3. Compare their outputs with the bundled ones

In the toy 7TM rhodopsin dataset, the two outgroup sequences are:

- `OG_WP_285271495`
- `OG_WP_136361479`

These IDs must be passed to the pipelines as the `--outgroup` argument.

### RAxML + PAML pipeline: quick start on the toy dataset

All commands below are meant to be run from `raxml_paml_pipeline/examples/`.
These commands are not meant to be an optimal analysis pipeline, but just to demonstrate how the scripts work on a small 7TM rhodopsin dataset.

**Run RAxML-NG ML tree inference** (using the settings used to generate the bundled toy example outputs):

```bash
raxml-ng \
  --all \
  --msa toy_7tm_rhodopsin.fasta \
  --model LG+F+G4 \
  --tree pars{10},rand{10} \
  --bs-trees 1000 \
  --bs-metric fbp,tbe \
  --seed 12345 \
  --outgroup OG_WP_285271495,OG_WP_136361479 \
  --prefix toy_7tm_rhodopsin
```
This produces, among others:

- `toy_7tm_rhodopsin.raxml.bestTree` – ML tree

**Run PAML codeml** for amino-acid ASR on the RAxML tree
(either via `run_paml_asr.py` or a manually edited codeml.ctl):

- `run_paml_asr.py`  
  Minimal helper script that generates a codeml `.ctl` file for an amino-acid LG+Γ ASR run and executes `codeml`. This is provided only as an **example** to reproduce the PAML settings used in the toy 7TM rhodopsin dataset. Most PAML users will prefer to edit `.ctl` files manually for their own data.
  
```bash
python run_paml_asr.py \
  --msa    toy_7tm_rhodopsin.fasta \
  --tree   toy_7tm_rhodopsin.raxml.bestTree \
  --prefix ASR_toy_7tm_rhodopsin
```
This produces, among others:

- `rst` – PAML ASR output (site-wise posteriors)

**Run the indel-aware ASR wrapper** to obtain gap-aware ancestors:

```bash 
bash /path/to/run_indel_aware_paml.sh \
  --msa    toy_7tm_rhodopsin.fasta \
  --tree   toy_7tm_rhodopsin.raxml.bestTree \
  --rst    rst \
  --prefix ASR_toy_7tm_rhodopsin \
  --outgroup "OG_WP_285271495,OG_WP_136361479"
```

Alternatively, copy the scripts in the `raxml_paml_pipeline/indel_aware` folder to the current working directory.

```bash 
bash run_indel_aware_paml.sh \
  --msa    toy_7tm_rhodopsin.fasta \
  --tree   toy_7tm_rhodopsin.raxml.bestTree \
  --rst    rst \
  --prefix ASR_toy_7tm_rhodopsin \
  --outgroup "OG_WP_285271495,OG_WP_136361479"
```
This produces, for all internal nodes:

- `ASR_toy_7tm_rhodopsin_indel_withgap.fasta` – indel-aware ancestor sequences retaining alignment gaps
- `ASR_toy_7tm_rhodopsin_indel_nogap.fasta` – gap-stripped versions used for AlphaFold and summary statistics
- `ASR_toy_7tm_rhodopsin_raw_asr_withgap.fasta` – raw ancestor sequences extracted directly from the PAML `rst` file retaining alignment gaps
- `ASR_toy_7tm_rhodopsin_raw_asr_nogap.fasta` – gap-stripped versions used for AlphaFold and summary statistics
- `ASR_toy_7tm_rhodopsin_paml_nodes.tree` – PAML-compatible tree with numeric internal-node labels; use this file to identify the PAML node IDs used in the output ancestral FASTA files
- `ASR_toy_7tm_rhodopsin_indel_work/` – intermediate files for binary gap-state reconstruction and node mapping

**Map PP / PP×pLDDT / PP−pLDDT onto an AlphaFold model** for a focal ancestor
(e.g. node name `28` in the toy example):

```bash
cd confmap/Node28   # directory containing the AlphaFold CIF/JSON
```
```bash
bash /path/to/run_confmap_paml.sh \
  --rst   /path/to/rst \
  --node  28 \
  --withgap  /path/to/ASR_toy_7tm_rhodopsin_indel_withgap.fasta \
  --outdir confmap
```
Alternatively, copy the scripts in the `raxml_paml_pipeline/confmap` folder, `rst` and `withgap.fasta` files to the current working directory.

```bash
bash run_confmap_paml.sh \
  --rst   rst \
  --node  28 \
  --withgap  ASR_toy_7tm_rhodopsin_indel_withgap.fasta \
  --outdir confmap
```
This creates in `confmap/`:

- `Node28_plddt_bfactor.pdb` – pLDDT in B-factor
- `Node28_pp_bfactor.pdb` – PP in B-factor
- `Node28_ppminusplddt_bfactor.pdb` – scaled PP−pLDDT in B-factor
- `Node28_ppxplddt_bfactor.pdb` – PP×pLDDT in B-factor
- `Node28_conf_CA_stats.log` – summary statistics for site-wise PP and CA-based structure metrics
- `Node28_aligned_withgap.fasta` - aligned indel-aware ancestral sequence used for coordinate mapping
- `Node28_pp_values.txt` - site-wise ancestral posterior probabilities, after taking the maximum over states per site for the node
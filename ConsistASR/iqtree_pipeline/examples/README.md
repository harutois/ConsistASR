## Toy example: 7TM rhodopsin dataset

This directory contains a small 7TM rhodopsin toy dataset and bundled example outputs for the IQ-TREE-based ConsistASR workflow.

* `toy_7tm_rhodopsin.fasta`

Example outputs (trees, ASR state files, indel-aware FASTA, confmap PDBs) are provided under the corresponding `examples/` subdirectories, so users can:

1. Inspect intermediate files
2. Re-run the scripts step by step
3. Compare their outputs with the bundled ones

In the toy 7TM rhodopsin dataset, the two outgroup sequences are:

- `OG_WP_285271495`
- `OG_WP_136361479`

These IDs must be passed to the pipelines as the `--outgroup` argument.

### IQ-TREE pipeline: quick start on the toy dataset

All commands below are meant to be run from `iqtree_pipeline/examples/`.
These commands are not meant to be an optimal analysis pipeline, but just to demonstrate how the scripts work on a small 7TM rhodopsin dataset.

### Run IQ-TREE ML tree inference

To avoid overwriting the bundled expected outputs, create a separate working directory:

```bash
mkdir -p reproduce_Tree
cd reproduce_Tree
```

```bash
iqtree3 \
  -s ../toy_7tm_rhodopsin.fasta \
  -m Q.pfam+I+G4 \
  -bb 1000 \
  -alrt 1000 \
  -o OG_WP_285271495,OG_WP_136361479 \
  -pre toy_7tm_rhodopsin \
  -seed 12345
```
This produces, among others:

- `toy_7tm_rhodopsin.treefile` – ML tree with SH-aLRT / UFBoot
- `toy_7tm_rhodopsin.contree` – consensus tree with UFBoot2

### Run IQ-TREE ASR
To avoid overwriting the bundled expected outputs, create a separate working directory:

```bash
mkdir -p reproduce_ASR
cd reproduce_ASR
```

```bash
iqtree3 \
  -s ../toy_7tm_rhodopsin.fasta \
  -m Q.pfam+I+G4 \
  -te ../Tree/toy_7tm_rhodopsin.treefile \
  -o OG_WP_285271495,OG_WP_136361479 \
  --ancestral \
  -pre ASR_toy_7tm_rhodopsin \
  -seed 12345
```
This produces, among others:

- `ASR_toy_7tm_rhodopsin.state` – IQ-TREE ASR output (site-wise posteriors)

### Run the indel-aware ASR wrapper

To avoid overwriting the bundled expected outputs:

```bash
mkdir -p reproduce_indel_aware
cd reproduce_indel_aware
```

```bash 
bash /path/to/run_indel_aware_iqtree.sh \
  --msa    ../toy_7tm_rhodopsin.fasta \
  --tree   ../ASR/ASR_toy_7tm_rhodopsin.treefile \
  --state  ../ASR/ASR_toy_7tm_rhodopsin.state \
  --prefix ASR_toy_7tm_rhodopsin \
  --outgroup "OG_WP_285271495,OG_WP_136361479"
```

This produces, for all reconstructed internal nodes:

- `ASR_toy_7tm_rhodopsin_raw_asr_withgap.fasta` – raw IQ-TREE ancestral sequences retaining the alignment coordinate system
- `ASR_toy_7tm_rhodopsin_raw_asr_nogap.fasta` – raw ancestral sequences with gap characters removed
- `ASR_toy_7tm_rhodopsin_indel_withgap.fasta` – indel-aware ancestral sequences retaining the alignment coordinate system
- `ASR_toy_7tm_rhodopsin_indel_nogap.fasta` – gap-stripped indel-aware ancestral sequences for downstream applications
- `ASR_toy_7tm_rhodopsin_binary_posterior.tsv` – posterior probabilities of ancestral residue absence and presence
- `ASR_toy_7tm_rhodopsin_site_posterior.tsv` – site-wise integration of IQ-TREE amino-acid reconstruction and binary gap-state reconstruction
- `ASR_toy_7tm_rhodopsin_indel_nogap_posterior.tsv` – posterior information corresponding directly to residues retained in the final gap-stripped ancestral sequences

Intermediate files used for the binary gap-state reconstruction are retained in:

- `ASR_toy_7tm_rhodopsin_indel_work/`

Representative deterministic output files, particularly the final FASTA files, can be compared directly with the bundled expected outputs using cmp.

```bash
cmp ASR_toy_7tm_rhodopsin_indel_nogap.fasta \
../indel_aware/ASR_toy_7tm_rhodopsin_indel_nogap.fasta
```

No output from `cmp` indicates that the generated file is identical to the bundled expected output.

Return to the example directory:

```bash
cd ..
```

### Map PP / PP×pLDDT / PP−pLDDT onto an AlphaFold model
The bundled example uses IQ-TREE `Node10`.

To avoid modifying the bundled AlphaFold example directory, create a separate working directory:

```bash
mkdir -p reproduce_confmap_Node10
cd reproduce_confmap_Node10
```

Copy the bundled AlphaFold Server output files used for this example:

```bash
cp ../confmap/Node10/toy_7tm_rhodopsin_node_10_model_0.cif .
cp ../confmap/Node10/toy_7tm_rhodopsin_node_10_summary_confidences_0.json .
```

Run the IQ-TREE ConfMap wrapper:

```bash
bash /path/to/run_confmap_paml.sh \
  --rst ../ASR/ASR_toy_7tm_rhodopsin.state \
  --node Node10 \
  --withgap ../indel_aware/ASR_toy_7tm_rhodopsin_indel_withgap.fasta \
  --outdir confmap
```

This creates in `confmap/`:

- `Node10_plddt_bfactor.pdb` – pLDDT in B-factor
- `Node10_pp_bfactor.pdb` – PP in B-factor
- `Node10_ppminusplddt_bfactor.pdb` – scaled PP−pLDDT in B-factor
- `Node10_ppxplddt_bfactor.pdb` – PP×pLDDT in B-factor
- `Node10_conf_CA_stats.log` – summary statistics for site-wise PP and CA-based structure metrics
- `Node10_aligned_withgap.fasta` - aligned indel-aware ancestral sequence used for coordinate mapping
- `Node10_pp_values.txt` - site-wise ancestral posterior probabilities, after taking the maximum over states per site for the node

These PP values are amino-acid ASR posterior probabilities from the IQ-TREE `.state` file and are distinct from the binary residue-presence/residue-absence posterior probabilities generated by the indel-aware workflow.
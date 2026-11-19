## Toy example: 7TM rhodopsin dataset

This directory contains a small 7TM rhodopsin toy dataset and bundled example outputs for the RAxML-NG + PAML-based ConsistASR workflow.

* `toy_7tm_rhodopsin.fasta`

Bundled outputs for tree inference, PAML amino-acid ASR, indel-aware processing, and confidence mapping are provided in the corresponding subdirectories so that users can:

1. Inspect intermediate files
2. Re-run the workflow step by step
3. Compare representative deterministic outputs with the bundled expected files

In the toy 7TM rhodopsin dataset, the two outgroup sequences are:

- `OG_WP_285271495`
- `OG_WP_136361479`

These IDs must be passed to the relevant commands as the `--outgroup` argument.

### RAxML-NG + PAML pipeline: quick start on the toy dataset

The commands below assume that the current directory is:

```text
raxml_paml_pipeline/examples/
```

They are intended to demonstrate and reproduce the ConsistASR toy workflow rather than to define an optimal phylogenetic-analysis strategy for all datasets.

### Run RAxML-NG ML tree inference

To avoid overwriting the bundled expected outputs, create a separate working directory:

```bash
mkdir -p reproduce_Tree
cd reproduce_Tree
```

Run RAxML-NG using the settings used to generate the bundled toy example:

```bash
raxml-ng \
  --all \
  --msa ../toy_7tm_rhodopsin.fasta \
  --model LG+F+G4 \
  --tree pars{10},rand{10} \
  --bs-trees 1000 \
  --bs-metric fbp,tbe \
  --seed 12345 \
  --outgroup OG_WP_285271495,OG_WP_136361479 \
  --prefix toy_7tm_rhodopsin
```

This produces, among other files:

- `toy_7tm_rhodopsin.raxml.bestTree` – maximum-likelihood tree

Return to the example directory before continuing:

```bash
cd ..
```

### Run PAML amino-acid ASR

The bundled `run_paml_asr.py` helper can be used to reproduce the PAML analysis used for the toy example. It prepares the required PAML input files and uses the bundled LG empirical amino-acid rate matrix by default.

Create a separate working directory:

```bash
mkdir -p reproduce_ASR
cd reproduce_ASR
```

Run:

```bash
python ../run_paml_asr.py \
  --msa ../toy_7tm_rhodopsin.fasta \
  --tree ../Tree/toy_7tm_rhodopsin.raxml.bestTree \
  --outdir paml_asr \
  --prefix ASR_toy_7tm_rhodopsin
```

The `--outdir` option keeps the reproduced PAML files separate from the bundled expected outputs.

Among the generated files:

- `paml_asr/rst` – PAML ancestral reconstruction output containing the node-labelled tree, reconstructed ancestral amino-acid states, and site-wise posterior probabilities

The helper also accepts `--aa-ratefile` when an alternative PAML-compatible empirical amino-acid rate matrix is required. If this option is omitted, the bundled LG rate matrix is used.

Return to the example directory:

```bash
cd ..
```

### Run the indel-aware ASR wrapper

To avoid overwriting the bundled expected outputs:

```bash
mkdir -p reproduce_indel_aware
cd reproduce_indel_aware
```

The following example uses the bundled PAML and tree outputs as fixed upstream inputs so that the indel-aware processing step can be reproduced independently:

```bash
bash /path/to/run_indel_aware_paml.sh \
  --msa ../toy_7tm_rhodopsin.fasta \
  --tree ../Tree/toy_7tm_rhodopsin.raxml.bestTree \
  --rst ../ASR/paml_asr/rst \
  --prefix ASR_toy_7tm_rhodopsin \
  --outgroup "OG_WP_285271495,OG_WP_136361479"
```

For this 18-sequence toy dataset, PAML reconstructs 17 ancestral nodes. The explicit PAML root is retained in the raw PAML outputs, whereas 16 non-root ancestral nodes are used for indel-aware integration.

Principal outputs are:

- `ASR_toy_7tm_rhodopsin_raw_asr_withgap.fasta` – raw PAML ancestral sequences retaining the alignment coordinate system
- `ASR_toy_7tm_rhodopsin_raw_asr_nogap.fasta` – raw PAML ancestral sequences with gap characters removed
- `ASR_toy_7tm_rhodopsin_indel_withgap.fasta` – indel-aware sequences for the 16 mapped non-root ancestral nodes, retaining alignment coordinates
- `ASR_toy_7tm_rhodopsin_indel_nogap.fasta` – final gap-stripped indel-aware ancestral sequences
- `ASR_toy_7tm_rhodopsin_binary_posterior.tsv` – posterior probabilities of ancestral residue absence and presence
- `ASR_toy_7tm_rhodopsin_site_posterior.tsv` – site-wise integration of PAML amino-acid reconstruction and binary gap-state information
- `ASR_toy_7tm_rhodopsin_indel_nogap_posterior.tsv` – posterior information corresponding directly to residues retained in the final gap-stripped ancestral sequences
- `ASR_toy_7tm_rhodopsin_paml_nodes.tree` – PAML node-labelled tree extracted from `rst`
- `ASR_toy_7tm_rhodopsin_node_map.tsv` – correspondence between the 16 non-root PAML ancestral nodes and RAxML-NG internal vertices

Intermediate files used for binary gap-state reconstruction and validation are retained in:

- `ASR_toy_7tm_rhodopsin_indel_work/`

Representative deterministic outputs can be compared directly with the bundled expected files. For example:

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

The bundled example uses PAML node `28`.

To avoid modifying the bundled AlphaFold example directory, create a separate working directory:

```bash
mkdir -p reproduce_confmap_Node28
cd reproduce_confmap_Node28
```

Copy the bundled AlphaFold Server output files used for this example:

```bash
cp ../confmap/Node28/toy_7tm_rhodopsin_rmp_node_28_model_0.cif .
cp ../confmap/Node28/toy_7tm_rhodopsin_rmp_node_28_summary_confidences_0.json .
```

Run the PAML ConfMap wrapper:

```bash
bash /path/to/run_confmap_paml.sh \
  --rst ../ASR/paml_asr/rst \
  --node 28 \
  --withgap ../indel_aware/ASR_toy_7tm_rhodopsin_indel_withgap.fasta \
  --outdir confmap
```

This creates in `confmap/`:

- `Node28_plddt_bfactor.pdb` – AlphaFold pLDDT stored in the B-factor field
- `Node28_pp_bfactor.pdb` – PAML amino-acid ASR posterior probability stored in the B-factor field
- `Node28_ppminusplddt_bfactor.pdb` – scaled PP − pLDDT
- `Node28_ppxplddt_bfactor.pdb` – PP × pLDDT
- `Node28_conf_CA_stats.log` – summary statistics for site-wise PP and CA-based structure metrics
- `Node28_aligned_withgap.fasta` – aligned ancestral sequence used for coordinate mapping
- `Node28_pp_values.txt` – site-wise posterior probabilities for PAML node 28 extracted from the `rst` output

For PAML, the posterior probability corresponding to the requested ancestral node is extracted directly from each site-wise record; no additional maximum-over-states aggregation is performed.

These PAML posterior probabilities describe confidence in the reconstructed amino-acid state and are distinct from the binary residue-presence/residue-absence posterior probabilities generated by the indel-aware workflow.
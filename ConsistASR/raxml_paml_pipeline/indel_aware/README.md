# Indel-aware ancestral sequence reconstruction (RAxML-NG + PAML)

This directory provides a small post-processing workflow that combines:

- **RAxML-NG** phylogeny used as the fixed topology for PAML amino-acid ASR
- **PAML** amino-acid ancestral sequence reconstruction (`rst`)
- **RAxML-NG** binary residue-presence/residue-absence reconstruction

to generate both **raw PAML ancestral sequences** and **indel-aware ancestral amino-acid sequences** in FASTA format.

The workflow produces:

- raw PAML ancestral sequences parsed directly from the `rst` file
- indel-aware ancestral sequences retaining alignment coordinates
- gap-stripped indel-aware ancestral sequences for downstream applications
- site-wise posterior probabilities for binary residue presence/absence
- integrated site-wise posterior information for the final ancestral sequences

This is the RAxML-NG + PAML indel-aware workflow of *ConsistASR*.

---

## Contents

Scripts in this directory:

- `run_indel_aware_paml.sh`  
  Main wrapper script that runs the complete RAxML-NG + PAML indel-aware post-processing workflow.

- `map_raxml_to_paml_nodes_from_rst.py`  
  Extracts the PAML node-labelled tree from `rst`, normalizes PAML-added numeric prefixes on tip labels when appropriate, and identifies correspondence between non-root PAML ancestral nodes and RAxML-NG internal vertices from their phylogenetic positions.

- `paml_state_and_indel_to_fasta.py`  
  Extracts raw PAML ancestral sequences and integrates the amino-acid reconstruction with the mapped binary residue-presence/residue-absence reconstruction.

Shared helper scripts are located in `../../common/`:

- `msa_to_binary.py`
- `raxmlng_binary_asr.py`

---

## Software requirements

- **Python 3.10+**
- **RAxML-NG**
  - used for binary-model evaluation and ancestral binary-state reconstruction
- **PAML**
  - PAML ASR is assumed to have been completed beforehand, producing an `rst` file containing the node-labelled tree and ancestral amino-acid states
- A POSIX-like shell environment such as Linux or macOS

> **Note:** PAML itself is not called by this wrapper. The workflow starts from an existing PAML `rst` file and the corresponding phylogeny used for that ASR analysis.

---

## What the pipeline does

### 1. Convert the amino-acid MSA to binary residue-presence/residue-absence states

The shared `common/msa_to_binary.py` script converts the aligned amino-acid sequences to a binary PHYLIP alignment.

`-` is encoded as residue absence (`0`). Standard and ambiguous amino-acid symbols representing an existing residue are encoded as presence (`1`). Unsupported gap-like or special characters are rejected.

### 2. Re-estimate branch lengths and the binary model on the fixed topology

RAxML-NG `--evaluate` is run with `BIN+G`.

Only the topology of the amino-acid phylogeny is retained. Branch lengths and the free parameters of the binary model are re-estimated from the binary residue-presence/residue-absence alignment.

Representative outputs include:

```text
<prefix>_indel_eval.raxml.bestModel
<prefix>_indel_eval.raxml.bestTree
```

### 3. Reconstruct ancestral binary states with RAxML-NG

RAxML-NG `--ancestral` is run using the evaluated tree and model, with model parameters and branch lengths fixed.

Representative outputs include:

```text
<prefix>_indel_asr.raxml.ancestralProbs
<prefix>_indel_asr.raxml.ancestralStates
<prefix>_indel_asr.raxml.ancestralTree
```

The shared `common/raxmlng_binary_asr.py` script validates the RAxML-NG probability and state files and derives binary residue-presence/residue-absence calls from the posterior probabilities.

The corresponding posterior table is written as:

```text
<prefix>_binary_posterior.tsv
```

### 4. Map non-root PAML ancestral nodes to RAxML-NG internal vertices

`map_raxml_to_paml_nodes_from_rst.py` extracts the PAML node-labelled tree from the `rst` file and compares phylogenetic positions between the PAML and RAxML-NG trees.

PAML represents a rooted tree and therefore contains one explicit degree-2 root in addition to the non-root internal nodes represented by the corresponding unrooted RAxML-NG topology.

The workflow therefore:

1. retains the explicit PAML root in the raw PAML ASR output;
2. excludes that root from binary gap-state integration;
3. requires all non-root PAML ancestral nodes to map one-to-one to RAxML-NG internal vertices.

For the bundled 18-sequence toy dataset, PAML reconstructs 17 ancestral nodes, whereas 16 non-root PAML nodes correspond to the 16 RAxML-NG internal vertices.

The mapping is written to:

```text
<prefix>_node_map.tsv
```

The extracted PAML node-labelled tree is written to:

```text
<prefix>_paml_nodes.tree
```

### 5. Extract raw PAML sequences and generate indel-aware ancestral sequences

Raw PAML ancestral sequences are extracted directly from `rst`:

```text
<prefix>_raw_asr_withgap.fasta
<prefix>_raw_asr_nogap.fasta
```

These raw files retain all ancestral nodes reconstructed by PAML, including its explicit root.

Binary gap-state information is then integrated for the mapped non-root nodes to produce:

```text
<prefix>_indel_withgap.fasta
<prefix>_indel_nogap.fasta
```

The explicit PAML root is not included in these indel-aware FASTA files because there is no one-to-one corresponding RAxML-NG internal vertex for binary-state integration.

Additional posterior tables are written as:

```text
<prefix>_site_posterior.tsv
<prefix>_indel_nogap_posterior.tsv
```

`*_site_posterior.tsv` provides a complete alignment-site audit trail. The explicit PAML root is retained in this table with binary-state fields unavailable because binary integration is not performed for that node.

`*_indel_nogap_posterior.tsv` contains posterior information corresponding directly to residues retained in the final gap-stripped indel-aware ancestral sequences and therefore contains only the mapped non-root nodes.

---

## Output summary

For the bundled 18-sequence toy dataset:

- `*_raw_asr_*.fasta` contains **17 PAML ancestral sequences**, including the explicit PAML root.
- `*_indel_*.fasta` contains **16 indel-aware ancestral sequences**, corresponding to the mapped non-root nodes.

Principal outputs are:

```text
<prefix>_raw_asr_withgap.fasta
<prefix>_raw_asr_nogap.fasta
<prefix>_indel_withgap.fasta
<prefix>_indel_nogap.fasta
<prefix>_binary_posterior.tsv
<prefix>_site_posterior.tsv
<prefix>_indel_nogap_posterior.tsv
<prefix>_paml_nodes.tree
<prefix>_node_map.tsv
```

Intermediate RAxML-NG and LF-normalized working files are retained in:

```text
<prefix>_indel_work/
```

---

## Example folder structure

A typical run may contain:

```text
.
├── HeR_SzR_228_FFT_NS_1.fasta
├── HeR_SzR_228_FFT_NS_1.raxml.bestTree
├── HeR_SzR_228_FFT_NS_1.rst
├── ASR_LGFG4_FFT_NS_1_raw_asr_withgap.fasta
├── ASR_LGFG4_FFT_NS_1_raw_asr_nogap.fasta
├── ASR_LGFG4_FFT_NS_1_indel_withgap.fasta
├── ASR_LGFG4_FFT_NS_1_indel_nogap.fasta
├── ASR_LGFG4_FFT_NS_1_binary_posterior.tsv
├── ASR_LGFG4_FFT_NS_1_site_posterior.tsv
├── ASR_LGFG4_FFT_NS_1_indel_nogap_posterior.tsv
├── ASR_LGFG4_FFT_NS_1_paml_nodes.tree
├── ASR_LGFG4_FFT_NS_1_node_map.tsv
└── ASR_LGFG4_FFT_NS_1_indel_work/
    ├── input_msa.lf.fasta
    ├── input_tree.lf.nwk
    ├── input_paml_rst.lf
    ├── ASR_LGFG4_FFT_NS_1_binary.phy
    ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.bestModel
    ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.bestTree
    ├── ASR_LGFG4_FFT_NS_1_indel_asr.raxml.ancestralProbs
    ├── ASR_LGFG4_FFT_NS_1_indel_asr.raxml.ancestralStates
    └── ASR_LGFG4_FFT_NS_1_indel_asr.raxml.ancestralTree
```

RAxML-NG may generate additional log or intermediate files depending on version and run conditions.

---

## Log messages and validation

A successful run should report that:

- the RAxML-NG binary probability and state files are internally consistent;
- all non-root PAML ancestral nodes map one-to-one to RAxML-NG internal vertices;
- the explicit PAML root is identified separately and excluded from indel-aware integration;
- the final indel-aware sequences and posterior tables are written successfully.

For the bundled toy dataset, 17 PAML ancestral nodes are expected in the raw PAML reconstruction, of which 16 non-root nodes are used for indel-aware integration.

An unexpected failure to map any non-root PAML ancestral node is treated as an error rather than being silently ignored.

---

## Limitations and scope

This workflow combines:

- PAML amino-acid ancestral reconstruction from `rst`
- RAxML-NG binary residue-presence/residue-absence reconstruction

The workflow assumes that the MSA, PAML `rst`, and RAxML-NG tree describe the same taxon set and underlying topology.

Because PAML represents an explicit rooted tree whereas RAxML-NG represents the corresponding unrooted internal vertices, the explicit PAML root is not assigned a binary gap-state reconstruction in the current workflow.

Binary gap-state reconstruction can be sensitive to evolutionary-model assumptions, particularly for deep ancestral nodes and gap-rich regions. Sensitivity analyses using alternative reasonable binary models should be considered for biologically important deep reconstructions.

The included examples focus on 7TM microbial rhodopsins, but the workflow can be applied to other sufficiently alignable protein families containing lineage-dependent insertions and deletions.

---

## Citation

If you use this workflow in a publication, please cite:

- PAML
- RAxML-NG
- the associated *ConsistASR* publication

The repository-level README provides the corresponding references and software archive information.
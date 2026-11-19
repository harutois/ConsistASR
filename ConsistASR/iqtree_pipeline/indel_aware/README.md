# Indel-aware ancestral sequence reconstruction (IQ-TREE + RAxML-NG)

This directory provides a small post-processing workflow that combines:

- **IQ-TREE** amino-acid ancestral sequence reconstruction (`.state`, `.treefile`)
- **RAxML-NG** binary residue-presence/residue-absence reconstruction

to generate both **raw IQ-TREE ancestral sequences** and **indel-aware ancestral amino-acid sequences** in FASTA format.

The workflow produces:

- raw IQ-TREE ancestral sequences before indel-aware processing
- indel-aware ancestral sequences retaining alignment coordinates
- gap-stripped indel-aware ancestral sequences for downstream applications
- site-wise posterior probabilities for binary residue presence/absence
- integrated site-wise posterior information for the final ancestral sequences

This is the IQ-TREE-based indel-aware workflow of *ConsistASR*.

---

## Contents

Scripts in this directory:

- `run_indel_aware_iqtree.sh`  
  Main wrapper script that runs the complete IQ-TREE indel-aware post-processing workflow.

- `state_and_indel_to_fasta.py`  
  Extracts raw IQ-TREE ancestral FASTA sequences and integrates the amino-acid reconstruction with the binary gap-state reconstruction.

Shared helper scripts are located in `../../common/`:

- `msa_to_binary.py`  
  Converts an amino-acid MSA into a binary residue-presence/residue-absence alignment for RAxML-NG.

- `raxmlng_binary_asr.py`  
  Parses and validates RAxML-NG binary ancestral-state probabilities and generates binary profiles and site-wise posterior output.

Within the standard repository layout, the wrapper automatically uses the shared helpers in `common/`. Copies placed alongside the wrapper can also be used for standalone distributions.

---

## Software requirements

- **Python 3.10+**
- **RAxML-NG**
  - Used for both binary-model evaluation and ancestral binary-state reconstruction.
- A POSIX-like shell environment (Linux / macOS)
  - `bash`, `mkdir`, `cp`, etc.

The reference ConsistASR conda environment is provided in the top-level `environment.yml`.

> **Note:** IQ-TREE itself is not called by this wrapper. The workflow assumes that IQ-TREE amino-acid ASR has already been completed and that the corresponding `.treefile` and `.state` files are available.

---

## Input files

The wrapper expects the IQ-TREE ASR outputs together with the exact amino-acid MSA used for that analysis.

### 1. Amino-acid MSA in FASTA format

The MSA must be the same alignment used for the IQ-TREE ASR analysis.

Do not realign, trim, reorder alignment columns, or otherwise modify the alignment before running the indel-aware workflow.

Example:

```text
HeR_SzR_228_PSI_TM_Coffee.fasta
```

For binary conversion:

- `-` is interpreted as residue absence (`0`)
- the 20 standard amino acids are interpreted as residue presence (`1`)
- `B`, `Z`, `J`, and `X` are also treated as residue presence because they indicate an existing residue with ambiguous amino-acid identity
- lowercase amino-acid letters are accepted
- other gap-like or special characters such as `.`, `?`, `~`, `*`, and `!` are rejected with an error

### 2. IQ-TREE tree file (`.treefile`)

The tree used for IQ-TREE ASR, containing the internal `NodeXX` labels associated with the `.state` file.

Example:

```text
ASR_QpfamR7_PSI_TM_Coffee.treefile
```

### 3. IQ-TREE ancestral-state file (`.state`)

The IQ-TREE site-wise ancestral amino-acid reconstruction.

Example:

```text
ASR_QpfamR7_PSI_TM_Coffee.state
```

### 4. Outgroup list

Comma-separated list of one or more taxa:

```text
OG_WP_285271495,OG_WP_136361479,OG_WP_010903286
```

Names must exactly match sequence IDs in the MSA.

---

## Usage

Basic example:

```bash
bash /path/to/run_indel_aware_iqtree.sh \
  --msa HeR_SzR_228_PSI_TM_Coffee.fasta \
  --tree ASR_QpfamR7_PSI_TM_Coffee.treefile \
  --state ASR_QpfamR7_PSI_TM_Coffee.state \
  --prefix ASR_QpfamR7_PSI_TM_Coffee \
  --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
```

### Command-line options

- `--msa`  
  Input amino-acid MSA in FASTA format.

- `--tree`  
  IQ-TREE `.treefile` used for ASR.

- `--state`  
  IQ-TREE `.state` file.

- `--prefix`  
  Prefix for newly generated files.

- `--outgroup`  
  Comma-separated list of outgroup taxa. Names must exactly match sequence IDs in the MSA.

---

## What the pipeline does

`run_indel_aware_iqtree.sh` performs the following steps.

### 1. Normalize and validate the input files

LF-normalized working copies of the MSA, tree, and IQ-TREE `.state` file are created inside:

```text
<prefix>_indel_work/
```

The original input files are not modified.

The wrapper checks input consistency before producing final ancestral sequences.

### 2. Convert the amino-acid MSA to binary residue-presence/residue-absence states

Shared script:

```text
common/msa_to_binary.py
```

The aligned amino-acid sequences are converted to a binary PHYLIP alignment:

```text
-       -> 0
residue -> 1
```

Output:

```text
<prefix>_binary.phy
```

### 3. Re-estimate branch lengths and the binary model on the fixed topology

RAxML-NG is run in `--evaluate` mode using the binary alignment and the IQ-TREE topology:

```bash
raxml-ng --evaluate \
  --msa <prefix>_binary.phy \
  --msa-format PHYLIP \
  --model BIN+G \
  --tree <treefile> \
  --prefix <prefix>_indel_eval \
  --outgroup <OUTGROUPS> \
  --seed 12345
```

Only the topology of the amino-acid phylogeny is retained. Branch lengths and the free parameters of the binary model are re-estimated from the residue-presence/residue-absence data.

Representative outputs include:

```text
<prefix>_indel_eval.raxml.bestModel
<prefix>_indel_eval.raxml.bestTree
```

### 4. Reconstruct ancestral binary gap states with RAxML-NG

The evaluated binary model and tree are used for marginal ancestral-state reconstruction with model parameters and branch lengths fixed:

```bash
raxml-ng --ancestral \
  --msa <prefix>_binary.phy \
  --msa-format PHYLIP \
  --model <prefix>_indel_eval.raxml.bestModel \
  --tree <prefix>_indel_eval.raxml.bestTree \
  --prefix <prefix>_indel_asr \
  --outgroup <OUTGROUPS> \
  --opt-model off \
  --opt-branches off \
  --seed 12345
```

Representative outputs include:

```text
<prefix>_indel_asr.raxml.ancestralProbs
<prefix>_indel_asr.raxml.ancestralStates
<prefix>_indel_asr.raxml.ancestralTree
```

### 5. Validate the binary reconstruction

Shared script:

```text
common/raxmlng_binary_asr.py
```

The script validates:

- agreement between RAxML-NG `ancestralProbs` and `ancestralStates`
- node-label correspondence
- alignment-site coordinates
- continuity of site numbering
- posterior-probability ranges
- agreement between reported binary states and maximum-posterior states

For the IQ-TREE workflow, the RAxML-NG node labels and site coordinates are additionally required to match the IQ-TREE `.state` file exactly.

For each site:

```text
p_present > p_absent  -> present
p_absent > p_present  -> absent
```

An exact posterior tie is reported as ambiguous and conservatively retained as present for sequence integration. The posterior probabilities remain available for inspection.

The resulting posterior table is written as:

```text
<prefix>_binary_posterior.tsv
```

### 6. Extract raw IQ-TREE ASR sequences and integrate binary gap states

`state_and_indel_to_fasta.py` extracts the raw IQ-TREE amino-acid reconstruction and combines it with the binary residue-presence/residue-absence calls.

Raw ancestral sequences:

```text
<prefix>_raw_asr_withgap.fasta
<prefix>_raw_asr_nogap.fasta
```

Indel-aware ancestral sequences:

```text
<prefix>_indel_withgap.fasta
<prefix>_indel_nogap.fasta
```

For each alignment position:

- inferred binary absence -> final gap
- inferred binary presence -> retain the reconstructed amino-acid state
- if IQ-TREE reports `State="-"` but the binary reconstruction indicates presence, the highest-posterior amino acid in the IQ-TREE `.state` probability columns is used

The raw FASTA output preserves the original IQ-TREE reconstructed state and is not altered by the binary reconstruction.

Additional posterior tables are written as:

```text
<prefix>_site_posterior.tsv
<prefix>_indel_nogap_posterior.tsv
```

`*_site_posterior.tsv` provides an alignment-coordinate audit trail containing the amino-acid reconstruction, binary posterior probabilities, and final integrated state.

`*_indel_nogap_posterior.tsv` contains posterior information corresponding directly to residues retained in the final gap-stripped ancestral sequences.

---

## Output summary

After successful execution, the principal outputs are:

### FASTA files

- `<prefix>_raw_asr_withgap.fasta`  
  Raw IQ-TREE ancestral amino-acid reconstruction retaining alignment coordinates.

- `<prefix>_raw_asr_nogap.fasta`  
  Raw IQ-TREE ancestral reconstruction with gap characters removed.

- `<prefix>_indel_withgap.fasta`  
  Indel-aware ancestral sequences retaining the original alignment coordinate system.

- `<prefix>_indel_nogap.fasta`  
  Gap-stripped indel-aware ancestral protein sequences suitable for structure prediction, construct design, and other downstream analyses.

### Posterior-probability tables

- `<prefix>_binary_posterior.tsv`  
  Posterior probabilities of ancestral residue absence and presence.

- `<prefix>_site_posterior.tsv`  
  Site-wise integration of IQ-TREE amino-acid reconstruction and binary gap-state reconstruction.

- `<prefix>_indel_nogap_posterior.tsv`  
  Posterior information corresponding directly to residues retained in the final ungapped ancestral sequence.

### Intermediate working directory

```text
<prefix>_indel_work/
```

This directory retains the binary alignment, LF-normalized input copies, RAxML-NG evaluation output, and RAxML-NG ancestral-state output so that the binary reconstruction can be inspected independently.

---

## Example folder structure

A typical run may produce a directory like:

```text
.
├── ASR_QpfamR7_PSI_TM_Coffee.state
├── ASR_QpfamR7_PSI_TM_Coffee.treefile
├── HeR_SzR_228_PSI_TM_Coffee.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_raw_asr_withgap.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_raw_asr_nogap.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_indel_withgap.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_indel_nogap.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_binary_posterior.tsv
├── ASR_QpfamR7_PSI_TM_Coffee_site_posterior.tsv
├── ASR_QpfamR7_PSI_TM_Coffee_indel_nogap_posterior.tsv
└── ASR_QpfamR7_PSI_TM_Coffee_indel_work
    ├── input_msa.lf.fasta
    ├── input_tree.lf.nwk
    ├── input_iqtree_state.lf.state
    ├── ASR_QpfamR7_PSI_TM_Coffee_binary.phy
    ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.bestModel
    ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.bestTree
    ├── ASR_QpfamR7_PSI_TM_Coffee_indel_asr.raxml.ancestralProbs
    ├── ASR_QpfamR7_PSI_TM_Coffee_indel_asr.raxml.ancestralStates
    └── ASR_QpfamR7_PSI_TM_Coffee_indel_asr.raxml.ancestralTree
```

RAxML-NG may generate additional log or intermediate files depending on the version and run conditions. The listing above shows the principal files relevant to inspection of the ConsistASR workflow.

---

## Log messages and validation

A successful run should include messages indicating that:

```text
[INFO] RAxML-NG ancestralProbs / ancestralStates node labels match exactly
[INFO] IQ-TREE / RAxML-NG node labels and site coordinates match exactly
[INFO] Binary posterior rows validated
[DONE] Binary posterior TSV written
```

Exact wording and counts depend on the input dataset.

Unlike earlier ConsistASR versions, the current IQ-TREE workflow does not perform a separate RAxML-to-IQ-TREE internal-node mapping step. The fixed topology retains the IQ-TREE `NodeXX` labels, and the workflow requires exact node-label and site-coordinate correspondence between the RAxML-NG binary reconstruction and the IQ-TREE `.state` file.

A mismatch is treated as an error rather than being silently ignored.

If the workflow stops, inspect the reported node, site, or taxon mismatch before proceeding.

---

## Naming and formatting caveats

### Sequence IDs

Sequence IDs should be:

- unique
- free of whitespace
- preferably composed of letters, numbers, and underscores

Avoid spaces, pipes, and unusual punctuation when possible.

The `fasta_id_tools/` directory in the ConsistASR repository provides utilities for cleaning or shortening sequence identifiers.

### Outgroup names

Outgroup names must exactly match sequence IDs in the MSA.

If an outgroup taxon is not found, the wrapper terminates with an error.

### Alignment characters

Only `-` is interpreted as residue absence.

Ambiguous amino-acid symbols such as `B`, `Z`, `J`, and `X` indicate residue presence and are therefore encoded as `1`.

Other gap-like or special characters are rejected rather than being silently interpreted.

---

## Limitations and scope

This workflow performs indel-aware post-processing by combining:

- IQ-TREE amino-acid ancestral reconstruction
- RAxML-NG binary residue-presence/residue-absence reconstruction

The binary analysis retains the amino-acid tree topology but re-estimates branch lengths and binary-model parameters from the binary alignment.

The workflow assumes that:

- the MSA, IQ-TREE tree, and IQ-TREE `.state` file originate from the same analysis
- the taxon sets are identical
- the internal `NodeXX` labels and alignment-site coordinates are consistent

Binary gap-state reconstruction can be sensitive to evolutionary-model assumptions, particularly for deep ancestral nodes and gap-rich regions. For biologically important deep reconstructions, sensitivity to alternative reasonable binary models should be considered.

ConsistASR does not resolve uncertainty in the underlying multiple-sequence alignment itself. Highly divergent or poorly alignable protein families therefore require separate assessment of alignment uncertainty.

The included examples focus on 7TM microbial rhodopsins, but the workflow can be applied to other sufficiently alignable protein families containing lineage-dependent insertions and deletions.

---

## Troubleshooting

### Line endings

The wrapper internally creates LF-normalized copies of the input MSA, tree, and IQ-TREE `.state` file inside:

```text
<prefix>_indel_work/
```

The original input files remain unchanged.

This prevents hidden Windows-style carriage returns from causing taxon-name or parsing problems.

### Outgroup taxon not found

If the wrapper reports that an outgroup taxon is absent, confirm the identifiers in the MSA, for example:

```bash
seqkit seq -n input.fasta
```

and compare them with the names supplied to `--outgroup`.

### Node or site mismatch

If the wrapper reports an IQ-TREE / RAxML-NG node-label or site-coordinate mismatch, verify that the MSA, `.treefile`, and `.state` file were generated from the same IQ-TREE analysis and were not subsequently modified.

The workflow intentionally stops rather than attempting to infer correspondence from incompatible inputs.

---

## Citation

If you use this workflow in a publication, please cite:

- IQ-TREE
- RAxML-NG
- the associated *ConsistASR* publication

The repository-level README provides the corresponding references and software archive information.
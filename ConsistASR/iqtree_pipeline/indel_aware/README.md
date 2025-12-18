# Indel-aware ancestral sequence reconstruction (IQ-TREE + RAxML)

This directory provides a small pipeline that combines:

- **IQ-TREE** amino-acid ancestral state reconstruction (`.state`, `.treefile`)
- **RAxML(-NG + RAxML-HPC)** binary (0/1) indel reconstruction

to generate **indel-aware ancestral amino-acid sequences** in FASTA format:

- an alignment-like version **with gaps** (for inspection)
- a **gap-stripped** version (e.g. for AlphaFold input)

This is the IQ-TREE-based part of the *ConsistASR* toolkit.

---

## Contents

Scripts in this directory:

- `run_indel_aware_iqtree.sh`  
  Main wrapper script that runs the full pipeline end-to-end.

- `msa_to_binary.py`  
  Converts an amino-acid MSA in FASTA format to a binary (0/1) MSA in PHYLIP format for RAxML.

- `map_raxml_to_iqtree_nodes.py`  
  Maps RAxML internal node labels (numeric) to IQ-TREE `NodeXX` labels by comparing **leaf sets** (clades).

- `state_and_indel_to_fasta.py`  
  Merges IQ-TREE `.state` residues and RAxML indel (0/1) patterns into indel-aware ancestral sequences and writes FASTA files.

---

## Software requirements

- **Python 3.7+**
  - Only standard library is used in this subpipeline (no extra packages required).

- **RAxML-NG**  
  - Used in `--evaluate` mode on a binary (0/1) alignment (`BIN+G`).

- **RAxML-HPC (classic RAxML)**  
  - Used with `-f A` and `-m BINGAMMA` to perform binary (indel) ASR.

- A POSIX-like shell environment (Linux / macOS)  
  - `bash`, `mkdir`, `mv` etc.

> **Note:** IQ-TREE itself is not called here; this pipeline assumes that IQ-TREE has already been run and produced the `.treefile` and `.state` files.

---

## Input files

The main wrapper expects IQ-TREE ASR outputs plus the original MSA:

1. **MSA in FASTA format**

   - Same alignment used for IQ-TREE ASR  
   - Example:  
     `HeR_SzR_228_PSI_TM_Coffee.fasta`

2. **IQ-TREE tree file (`.treefile`)**

   - Tree used for ASR (typically containing `NodeXX` internal labels).  
   - Example:  
     `ASR_QpfamR7_PSI_TM_Coffee.treefile`

3. **IQ-TREE ancestral state file (`.state`)**

   - IQ-TREE amino-acid ancestral states at internal nodes.  
   - Example:  
     `ASR_QpfamR7_PSI_TM_Coffee.state`

4. **Outgroup list**

   - Comma-separated list of one or more taxa used as outgroup for RAxML-NG:
     ```
     OG_WP_285271495,OG_WP_136361479,OG_WP_010903286
     ```
   - **Names must exactly match** the sequence IDs in the MSA (after any cleaning/renaming you apply).

---

## Usage

Basic example:

```bash
bash run_indel_aware_iqtree.sh \
  --msa HeR_SzR_228_PSI_TM_Coffee.fasta \
  --tree ASR_QpfamR7_PSI_TM_Coffee.treefile \
  --state ASR_QpfamR7_PSI_TM_Coffee.state \
  --prefix ASR_QpfamR7_PSI_TM_Coffee \
  --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
```

### Command line options

* `--msa` (required)
  Input amino-acid MSA in FASTA format.

* `--tree` (required)
  IQ-TREE `.treefile` used for ASR.

* `--state` (required)
  IQ-TREE `.state` file.

* `--prefix` (required)
  Prefix for all newly generated files.

* `--outgroup` (required)
  Comma-separated list of outgroup taxa; must match sequence names in the MSA.

* `--workdir` (optional)
  Subdirectory name to store intermediate files.
  Default: `<prefix>_indel_work`

---

## What the pipeline does

`run_indel_aware_iqtree.sh` executes the following steps:

1. **MSA → binary 0/1 alignment (PHYLIP)**

   * Script: `msa_to_binary.py`
   * Rule: `'-' → 0`, any residue (`A–Z`) → `1`
   * Output:

     * `<prefix>_binary.phy` (in `<prefix>_indel_work/`)

2. **RAxML-NG likelihood evaluation on the binary data**

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

   * Output (in workdir): best-scoring tree

     * `<prefix>_indel_eval.raxml.bestTree`

3. **RAxML-HPC binary ASR (indel reconstruction)**

   ```bash
   raxmlHPC -f A \
     -m BINGAMMA \
     -t <prefix>_indel_eval.raxml.bestTree \
     -s <prefix>_binary.phy \
     -n <prefix>_binary \
     -F
   ```

   * Key outputs (in workdir):

     * `RAxML_nodeLabelledRootedTree.<prefix>_binary`
     * `RAxML_marginalAncestralStates.<prefix>_binary`

4. **Map RAxML node IDs to IQ-TREE `NodeXX`**

   * Script: `map_raxml_to_iqtree_nodes.py`
   * Compares **leaf sets** of each internal node between:

     * IQ-TREE ASR tree (`.treefile`)
     * RAxML node-labelled tree (`RAxML_nodeLabelledRootedTree.*`)
   * Outputs:

     * Indel file with IQ-TREE-style node names:
       `<prefix>_indel_ASR.iqtree_named.txt`
     * Node mapping table (TSV):
       `<prefix>_node_map.tsv`
       with columns:

       ```text
       iqtree_node   raxml_node   n_tips   status
       ```

5. **Merge IQ-TREE `.state` and RAxML indel data**

   * Script: `state_and_indel_to_fasta.py`
   * Rule for each position:

     * `bit = 0` → `'-'` (deletion)
     * `bit = 1` → keep amino-acid state (if AA is already `'-'`, keep `'-'`)
   * Outputs (in the **current directory**):

     * `<prefix>_indel_withgap.fasta`
       (alignment-style ancestral sequences, gaps preserved)
     * `<prefix>_indel_nogap.fasta`
       (all `'-'` removed, suitable as AlphaFold input)

6. **Move intermediate files**

   * All intermediate files are moved to `<prefix>_indel_work/`
   * Final outputs remain in the top-level directory:

     * `<prefix>_indel_withgap.fasta`
     * `<prefix>_indel_nogap.fasta`
     * `<prefix>_node_map.tsv`

---

## Output summary

After successful execution, you should see:

* **Main outputs**

  * `<prefix>_indel_withgap.fasta`
    Indel-aware ancestral sequences, alignment length identical to the MSA.
  * `<prefix>_indel_nogap.fasta`
    Gap-stripped ancestral sequences (one sequence per node), same `NodeXX` names.

* **Mapping table**

  * `<prefix>_node_map.tsv`
    Internal node mapping between RAxML and IQ-TREE:

    * `status = mapped` → safe 1:1 mapping
    * `status = unmapped` → typically a single node near a polytomy / outgroup

* **Intermediate files**

  * `<prefix>_indel_work/`
    Contains:

    * `<prefix>_binary.phy`, reduced alignments, RAxML log files,
      `RAxML_marginalAncestralStates.*`, `RAxML_nodeLabelledRootedTree.*`, etc.

---

## Example folder structure

A typical run might produce a directory like:

```text
test/
├── ASR_QpfamR7_PSI_TM_Coffee.state
├── ASR_QpfamR7_PSI_TM_Coffee.treefile
├── ASR_QpfamR7_PSI_TM_Coffee_indel_ASR.iqtree_named.txt
├── ASR_QpfamR7_PSI_TM_Coffee_indel_nogap.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_indel_withgap.fasta
├── ASR_QpfamR7_PSI_TM_Coffee_indel_work
│   ├── ASR_QpfamR7_PSI_TM_Coffee_binary.phy
│   ├── ASR_QpfamR7_PSI_TM_Coffee_binary.phy.reduced
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.bestModel
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.bestTree
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.bestTreeCollapsed
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.log
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.rba
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.reduced.phy
│   ├── ASR_QpfamR7_PSI_TM_Coffee_indel_eval.raxml.startTree
│   ├── RAxML_info.ASR_QpfamR7_PSI_TM_Coffee_binary
│   ├── RAxML_marginalAncestralProbabilities.ASR_QpfamR7_PSI_TM_Coffee_binary
│   ├── RAxML_marginalAncestralStates.ASR_QpfamR7_PSI_TM_Coffee_binary
│   └── RAxML_nodeLabelledRootedTree.ASR_QpfamR7_PSI_TM_Coffee_binary
├── ASR_QpfamR7_PSI_TM_Coffee_node_map.tsv
├── HeR_SzR_228_PSI_TM_Coffee.fasta
└── (pipeline scripts)
```

---

## Log messages and common warnings

Typical informational output:

```text
[INFO] Binary MSA written to: ASR_QpfamR7_PSI_TM_Coffee_indel_work/ASR_QpfamR7_PSI_TM_Coffee_binary.phy
...
[INFO] IQ-TREE internal nodes with NodeXX labels: 226
[INFO] RAxML internal nodes with numeric labels: 227
[INFO] Mapped internal nodes: 226
[WARN] 1 RAxML nodes could not be mapped to any NodeXX (topology mismatch or polytomy near outgroup?)
...
[WARN] 1 nodes found in .state but missing from indel file (skipped).
[WARN] 1 nodes found in indel file but not in .state (ignored).
[INFO] Merged 226 sequences
```

Notes:

* A small number of **unmapped nodes** (often exactly 1) is common when:

  * IQ-TREE and RAxML treat outgroup / polytomies slightly differently.
* These extra nodes are **excluded** from the final FASTA output.
* The number of merged sequences should match the number of **IQ-TREE internal nodes with `NodeXX` labels**.

If something looks suspicious:

* Inspect `<prefix>_node_map.tsv` and check which node(s) are `status = unmapped`.
* Confirm that all `NodeXX` of interest are present in `<prefix>_indel_nogap.fasta`.

---

## Naming and formatting caveats

* **Sequence IDs in the MSA**

  * PHYLIP format used by RAxML-NG/RAxML is strict:

    * No spaces in sequence names
    * 10-character limit in “classic” PHYLIP; here we use a simple, more permissive writer, but **avoid spaces and exotic symbols**.
  * It is strongly recommended to clean/normalize sequence IDs beforehand
    (e.g. replace spaces and special characters with underscores).

* **Outgroup names**

  * Must exactly match the sequence IDs used in the MSA.
  * If any outgroup taxon is not found in the alignment, the pipeline will stop with an error.

---

## Limitations and scope

* This pipeline is designed for **indel-aware ASR** using:

  * IQ-TREE amino-acid states
  * RAxML binary indel states
* It assumes that:

  * The IQ-TREE ASR tree and RAxML trees are based on **the same taxon set**.
  * Topologies are mostly consistent (except possibly near outgroup/polytomies).

It has been tested primarily on **7TM microbial rhodopsins**, but the approach is general and should be applicable to other protein families as long as these assumptions hold.

---

## Citation

If you use this pipeline in a publication, please cite:

* IQ-TREE / IQ-TREE3
* RAxML-NG
* RAxML (classic RAxML)

and the *ConsistASR* methodology paper (once available).

```text
[To be updated with the final ConsistASR reference]
```

---
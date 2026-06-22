# Indel-aware ancestral sequence reconstruction (RAxML + PAML)

This directory provides a small pipeline that combines:

- **RAxML** amino-acid phylogeny (used as the fixed tree for PAML ASR)
- **PAML** amino-acid ancestral sequence reconstruction (`rst`)
- **RAxML(-NG + RAxML-HPC)** binary (0/1) indel reconstruction

to generate both **raw PAML ASR sequences** and **indel-aware ancestral amino-acid sequences** in FASTA format:

- raw PAML ASR sequences parsed directly from the `rst` file
- indel-aware alignment-like sequences **with gaps** (for inspection)
- indel-aware **gap-stripped** sequences (e.g. for AlphaFold input)

This is the RAxML+PAML-based part of the *ConsistASR* toolkit.

---

## Contents

Scripts in this directory:

- `run_indel_aware_paml.sh`  
  Main wrapper script that runs the full RAxML+PAML indel-aware pipeline end-to-end.

- `msa_to_binary.py`  
  Converts an amino-acid MSA in FASTA format to a binary (0/1) MSA in PHYLIP format for RAxML.

- `map_raxml_to_paml_nodes_from_rst.py`  
  Extracts the PAML tree with node labels from a `rst` file, applies a simple “tail-normalization” to leaf names, and maps RAxML internal node labels to PAML node IDs by comparing **leaf sets** (clades).  
  It then rewrites the node labels in `RAxML_marginalAncestralStates.*` to the corresponding PAML node IDs.

- `paml_state_and_indel_to_fasta.py`  
  Parses PAML `rst` AA states for **all nodes**, writes raw PAML ASR FASTA files, then merges the AA states with RAxML indel (0/1) patterns already renamed to PAML node IDs.

---

## Software requirements

- **Python 3.7+**
  - Only the Python standard library is used in this subpipeline (no extra packages required).

- **RAxML-NG**  
  - Used in `--evaluate` mode on a binary (0/1) alignment (`BIN+G`).

- **RAxML-HPC (classic RAxML)**  
  - Used with `-f A` and `-m BINGAMMA` to perform binary (indel) ASR.

- **PAML**  
  - `codeml` (or similar) is assumed to have been run beforehand to produce a `rst` file containing:
    - the tree with node labels (“tree with node labels for Rod Page's TreeView”),
    - and per-node amino-acid states.

- A POSIX-like shell environment (Linux / macOS)  
  - `bash`, `mkdir`, `mv`, etc.

> **Note:** PAML itself is **not** called by this pipeline.  
> It assumes that PAML ASR has already been performed and that the `rst` file and the RAxML tree used for that ASR are available.

---

## Input files

The main wrapper expects:

1. **MSA in FASTA format**

   - Same alignment used for the RAxML+PAML ASR run.  
   - Example:  
     `HeR_SzR_228_FFT_NS_1.fasta`

2. **RAxML tree file used for PAML ASR**

   - RAxML tree given to PAML as the fixed topology.  
   - Example:  
     `HeR_SzR_228_FFT_NS_1.raxml.bestTree.tre`

3. **PAML ASR result file (`rst`)**

   - Contains:
     - “tree with node labels for Rod Page's TreeView”
     - “node # <id>” blocks with AA states for internal nodes.  
   - Example:  
     `HeR_SzR_228_FFT_NS_1.rst`

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
bash /path/to/run_indel_aware_paml.sh \
  --msa   HeR_SzR_228_FFT_NS_1.fasta \
  --tree  HeR_SzR_228_FFT_NS_1.raxml.bestTree.tre \
  --rst   HeR_SzR_228_FFT_NS_1.rst \
  --prefix ASR_LGFG4_FFT_NS_1 \
  --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
```

### Command line options

* `--msa` (required)
  Input amino-acid MSA in FASTA format.

* `--tree` (required)
  RAxML tree used for PAML ASR (Newick).

* `--rst` (required)
  PAML `rst` file containing the tree with node labels and ancestral states.

* `--prefix` (required)
  Prefix for all newly generated files.

* `--outgroup` (required)
  Comma-separated list of outgroup taxa; must match sequence names in the MSA.

*(A separate `--workdir` option is not exposed; intermediate files are automatically placed in `<prefix>_indel_work/`.)*

---

## What the pipeline does

`run_indel_aware_paml.sh` executes the following steps:

1. **MSA → binary 0/1 alignment (PHYLIP)**

   * Script: `msa_to_binary.py`
   * Rule: `'-' → 0`, any residue (`A–Z`) → `1`
   * Output:

     * `<prefix>_binary.phy` (later moved into `<prefix>_indel_work/`)

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

   * Output (in workdir): best-scoring tree on the binary data

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
       (0/1 indel states per node)

4. **Map RAxML node IDs to PAML node IDs**

   * Script: `map_raxml_to_paml_nodes_from_rst.py`

   * Procedure:

     1. Extract the PAML tree with node labels from the `rst` file.
     2. Apply *tail-normalization* on PAML leaf names:

        * Example:

          * PAML leaf: `112_OG_WP_010903286`
          * RAxML leaf: `OG_WP_010903286`
          * Normalization: `112_OG_WP_010903286 → OG_WP_010903286`
     3. Compare **leaf sets** (clades) between:

        * PAML tree (with PAML node IDs)
        * `RAxML_nodeLabelledRootedTree.*` (with RAxML numeric node labels)
     4. Build a mapping: **RAxML node → PAML node**.
     5. Rewrite the node labels in `RAxML_marginalAncestralStates.*` accordingly.

   * Outputs:

     * Indel file with **PAML node IDs**:
       `<prefix>_indel_ASR.paml_named.txt`
       (used as `--indel` input for the next step)
     * Node mapping table (TSV):
       `<prefix>_node_map.tsv`
       with columns:

       ```text
       paml_node   raxml_node   n_tips   status
       ```

5. **Write raw PAML ASR FASTA and merge `rst` AA states with RAxML indel data**

   * Script: `paml_state_and_indel_to_fasta.py`

   * Raw PAML ASR outputs are written directly from the `rst` file before indel correction:

     * `<prefix>_raw_asr_withgap.fasta`  
       Raw PAML ancestral sequences at alignment length.

     * `<prefix>_raw_asr_nogap.fasta`  
       Raw PAML ancestral sequences with all `'-'` characters removed.

   * Indel-aware correction rule for each position:

     * `bit = 0` → `'-'` (deletion at that column)
     * `bit = 1` → keep the AA state (if AA is already `'-'`, keep `'-'`)

   * Indel-aware outputs:

     * `<prefix>_indel_withgap.fasta`  
       Alignment-style ancestral sequences after indel correction.

     * `<prefix>_indel_nogap.fasta`  
       Gap-stripped indel-aware ancestral sequences, suitable for AlphaFold input.

   * Node-labelled tree output:

     * `<prefix>_paml_nodes.tree`  
       PAML node-labelled tree extracted from the `rst` file.

6. **Move intermediate files**

   * Intermediate files are moved to `<prefix>_indel_work/`, including:

     * `<prefix>_binary.phy`
     * `<prefix>_indel_eval.raxml.*`
     * `RAxML_*.<prefix>_binary`
     * `<prefix>_indel_ASR.paml_named.txt`

   * Final user-facing outputs remain in the top-level directory:

     * `<prefix>_raw_asr_withgap.fasta`
     * `<prefix>_raw_asr_nogap.fasta`
     * `<prefix>_indel_withgap.fasta`
     * `<prefix>_indel_nogap.fasta`
     * `<prefix>_paml_nodes.tree`
     * `<prefix>_node_map.tsv`

---

## Example folder structure

A typical run might produce a directory like:

```text
test/
├── ASR_LGFG4_FFT_NS_1_raw_asr_withgap.fasta
├── ASR_LGFG4_FFT_NS_1_raw_asr_nogap.fasta
├── ASR_LGFG4_FFT_NS_1_indel_nogap.fasta
├── ASR_LGFG4_FFT_NS_1_indel_withgap.fasta
├── ASR_LGFG4_FFT_NS_1_paml_nodes.tree
├── ASR_LGFG4_FFT_NS_1_node_map.tsv
├── ASR_LGFG4_FFT_NS_1_indel_work/
│   ├── ASR_LGFG4_FFT_NS_1_indel_ASR.paml_named.txt
│   ├── ASR_LGFG4_FFT_NS_1_binary.phy
│   ├── ASR_LGFG4_FFT_NS_1_binary.phy.reduced
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.bestModel
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.bestTree
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.bestTreeCollapsed
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.log
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.rba
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.reduced.phy
│   ├── ASR_LGFG4_FFT_NS_1_indel_eval.raxml.startTree
│   ├── RAxML_info.ASR_LGFG4_FFT_NS_1_binary
│   ├── RAxML_marginalAncestralProbabilities.ASR_LGFG4_FFT_NS_1_binary
│   ├── RAxML_marginalAncestralStates.ASR_LGFG4_FFT_NS_1_binary
│   └── RAxML_nodeLabelledRootedTree.ASR_LGFG4_FFT_NS_1_binary
├── HeR_SzR_228_FFT_NS_1.fasta
├── HeR_SzR_228_FFT_NS_1.raxml.bestTree.tre
├── HeR_SzR_228_FFT_NS_1.rst
└── (pipeline scripts)
```

---

## Log messages and common warnings

Typical informational output:

```text
[INFO] PAML leaves (raw) example: [...]
[INFO] PAML leaves (tail) example: [...]
[INFO] RAxML leaves example:       [...]
[INFO] Leaf-name overlap size:     228
[INFO] PAML internal nodes with labels: 227
[INFO] RAxML internal nodes with labels: 227
[INFO] Mapped internal nodes: 227
[INFO] Writing node mapping table to: ASR_LGFG4_FFT_NS_1_node_map.tsv
[DONE] ancestralStates written to: ASR_LGFG4_FFT_NS_1_indel_ASR.paml_named.txt
[INFO] Lines processed: 227, renamed: 2
...
[INFO] Reading PAML rst: HeR_SzR_228_FFT_NS_1.rst
[INFO] Loaded 227 nodes from rst
[INFO] Reading indel file: ASR_LGFG4_FFT_NS_1_indel_ASR.paml_named.txt
[INFO] Loaded 227 indel profiles
[INFO] Merging AA states and indel patterns...
[INFO] Merged 227 sequences
```

Notes:

* A small number of **renamed** node labels (e.g. `ROOT → 229`) is expected:

  * Most RAxML node numbers already match PAML’s node IDs.
  * Often only the root (or a few nodes) require renaming.
* If some nodes are missing in either the `rst` or the indel file, you may see warnings like:

  * `nodes found in rst but missing from indel file (skipped)`
  * `nodes found in indel file but not in rst (ignored)`

To inspect the mapping, open `<prefix>_node_map.tsv` and check:

* `status = mapped` → clade leaf set matches exactly between RAxML and PAML.
* `status = unmapped` → topology mismatch or polytomy in that clade.

---

## Naming and formatting caveats

* **Sequence IDs in the MSA**

  * PHYLIP format used by RAxML-NG/RAxML is strict:

    * No spaces in sequence names.
    * Classic PHYLIP requires short IDs; here we use a simple writer, but **avoid spaces and exotic characters** anyway.

  * It is strongly recommended to clean/normalize sequence IDs beforehand
    (e.g. replace spaces and special characters with underscores).

* **Outgroup names**

  * Must exactly match the sequence IDs used in the MSA.
  * If any outgroup taxon is not found in the alignment, the pipeline will stop with an error.

* **PAML leaf naming**

  * PAML often prefixes leaf names with `<index>_`, e.g.:

    ```text
    112_OG_WP_010903286
    ```

  * `map_raxml_to_paml_nodes_from_rst.py` assumes this format and strips the numeric prefix during tail-normalization so that leaves match the MSA / RAxML names (e.g. `OG_WP_010903286`).

---

## Limitations and scope

* This pipeline is designed for **indel-aware ASR** using:

  * PAML amino-acid states from `rst`
  * RAxML binary indel states

* It assumes that:

  * The PAML `rst` tree and the RAxML trees are based on **the same taxon set**.
  * Topologies are mostly consistent (except possibly near outgroup/polytomies).

It has been tested primarily on **7TM microbial rhodopsins**, but the approach is general and should be applicable to other protein families as long as these assumptions hold.

---

## Citation

If you use this pipeline in a publication, please cite:

* RAxML-NG
* RAxML (classic RAxML)
* PAML

and the *ConsistASR* paper.

---
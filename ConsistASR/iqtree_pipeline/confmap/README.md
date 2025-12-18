# confmap (IQ-TREE + AlphaFold) pipeline

This directory contains a small pipeline to map **ancestral sequence reconstruction (ASR)** posterior probabilities (PP) onto **AlphaFold** models as B-factors, and to derive derived metrics such as **PP × pLDDT** and **PP − pLDDT** for structural visualization and analysis.

The workflow is designed for **IQ-TREE-based ASR on indel-aware alignments**, combined with **AlphaFold monomer** predictions for selected ancestral nodes (e.g., Anc-SzR, Anc-HeR, Node112, etc.).

---

## 1. Files

- `run_confmap_iqtree.sh`  
  Main wrapper script. For a given internal node:
  - finds the corresponding AlphaFold model (`*_summary_confidences_*.json`, `*_model_X.cif`)
  - converts mmCIF → PDB (keeping pLDDT in B-factor)
  - extracts site-wise PP from the IQ-TREE `.state` file
  - maps PP to B-factors (per residue) using `map_confidence_to_bfactor.py`
  - generates additional PDB files with **PP − pLDDT (scaled)** and **PP × pLDDT**
  - computes CA-only statistics and writes a log file

- `map_confidence_to_bfactor.py`  
  Python script (using `gemmi`) that:
  - reads an AlphaFold mmCIF file
  - reads an aligned ASR sequence (with gaps)
  - reads per-site PP values
  - maps PP and/or pLDDT to per-residue B-factors, handling gaps explicitly

---

## 2. Dependencies

### Runtime

- Bash (POSIX shell; tested with `bash`)
- Python 3.x
- Python libraries:
  - `gemmi` (mmCIF parsing and PDB output)
  - `biopython` (`Bio.PDB` module for simple PDB manipulation)
- One of the following mmCIF → PDB converters available on `$PATH`:
  - `phenix.cif_as_pdb`
  - `gemmi` (command-line `gemmi convert`)
  - `pdb_fromcif`

### Upstream (assumed already done)

- **IQ-TREE** ASR:
  - `.state` file with site-wise posterior probabilities (e.g., `ASR_QpfamR7_PSI_TM_Coffee.state`)
  - **indel-aware alignment** with gaps, where the internal node of interest appears as a FASTA entry  
    (e.g., `ASR_QpfamR7_PSI_TM_Coffee_indel_withgap.fasta`)
- **AlphaFold** run for the ancestral sequence of a specific node (Node112, Anc-SzR, etc.), producing a directory like:

  ```text
  Node112/
    fold_YYYY_MM_DD_HH_MM_..._summary_confidences_0.json
    fold_YYYY_MM_DD_HH_MM_..._summary_confidences_1.json
    ...
    fold_YYYY_MM_DD_HH_MM_..._model_0.cif
    fold_YYYY_MM_DD_HH_MM_..._model_1.cif
    ...
    msas/
    templates/
    terms_of_use.md
````

---

## 3. Directory layout and basic idea

The intended usage is:

1. AlphaFold is run for an ancestral node (e.g., Node112).
   This creates a directory such as `Node112/` containing several `*_summary_confidences_*.json` files and `*_model_X.cif` files.

2. `run_confmap_iqtree.sh` is executed **inside that AlphaFold output directory**, using:

   * an IQ-TREE `.state` file
   * an indel-aware alignment FASTA file
   * the node ID (e.g., `Node112`)

3. The script creates a subdirectory (by default `confmap/`) and writes all derived PDBs and a log file there.

Typical directory:

```text
Node112/
  fold_..._node112_summary_confidences_0.json
  fold_..._node112_model_0.cif
  ...
  msas/
  templates/
  terms_of_use.md
  confmap/
    Node112_plddt_bfactor.pdb
    Node112_pp_bfactor.pdb
    Node112_ppminusplddt_bfactor.pdb
    Node112_ppxplddt_bfactor.pdb
    Node112_aligned_indel.fasta
    Node112_pp_values.txt
    Node112_conf_CA_stats.log
```

---

## 4. Inputs

### 4.1 IQ-TREE ASR outputs

* `--state <STATE_FILE>`
  IQ-TREE `.state` file containing site-wise posterior probabilities, e.g.:

  ```bash
  ASR_QpfamR7_PSI_TM_Coffee.state
  ```

* `--aln <INDEL_WITHGAP_FASTA>`
  Indel-aware alignment (FASTA) that includes the node of interest as a header:

  ```bash
  ASR_QpfamR7_PSI_TM_Coffee_indel_withgap.fasta
  ```

  The script expects an exact header `>${NODE}` (e.g., `>Node112`).

### 4.2 AlphaFold outputs

From the AlphaFold run for the ancestral sequence of the node:

* `*_summary_confidences_*.json`
  (optionally containing the node id in the filename, e.g. `_node112_`)
* `*_model_X.cif` (corresponding mmCIF models)

The script automatically:

1. Collects all `*_summary_confidences_*.json`
2. Prefers those whose filename contains the node id (case-insensitive)
3. Among the remaining candidates, selects the one with the highest `ranking_score` (if present) or `plddt`
4. Extracts the model index `X` from the filename suffix `_summary_confidences_X.json`
5. Uses the corresponding `*_model_X.cif` as the best AlphaFold model

---

## 5. Usage

Run **inside** the AlphaFold directory for a given node:

```bash
cd Node112/

bash /path/to/run_confmap_iqtree.sh \
  --state /path/to/ASR_QpfamR7_PSI_TM_Coffee.state \
  --node  Node112 \
  --aln   /path/to/ASR_QpfamR7_PSI_TM_Coffee_indel_withgap.fasta \
  --outdir confmap
```

Arguments:

* `--state` : path to IQ-TREE `.state` file
* `--node`  : node ID (must match the FASTA header in `--aln`)
* `--aln`   : aligned FASTA (`with gaps`) containing that node sequence
* `--outdir`: output directory for confmap results (default: `confmap`)

Help:

```bash
bash run_confmap_iqtree.sh --help
```

---

## 6. Outputs

Within `--outdir` (default `confmap/`), the script writes:

* `NODE_plddt_bfactor.pdb`

  * PDB converted from the selected AlphaFold `.cif`
  * B-factors = original pLDDT values by AlphaFold

* `NODE_pp_bfactor.pdb`

  * PDB where B-factors = ASR posterior probability (PP) per residue
  * PP is automatically scaled:

    * if max(PP) ≤ 1.0, treat as 0–1 and scale to 0–100
    * otherwise assume already in 0–100

* `NODE_ppminusplddt_bfactor.pdb`

  * B-factor encodes **scaled (PP − pLDDT)**:

    * raw difference is `diff = PP − pLDDT` (range approx −100 to +100)
    * then linearly scaled to `[0, 100]` via `(diff + 100)/200 × 100`

* `NODE_ppxplddt_bfactor.pdb`

  * B-factor encodes **PP × pLDDT / 100**, i.e. a simple product metric:

    * PP and pLDDT are treated as 0–100 values
    * product is `(PP * pLDDT) / 100`, constrained to 0–100

* `NODE_aligned_indel.fasta`

  * the extracted aligned sequence for this node (one entry, with gaps)

* `NODE_pp_values.txt`

  * site-wise PP values (one per alignment column) extracted from `.state`, **after taking the maximum over states per site** for the node

* `NODE_conf_CA_stats.log`

  * CA-only statistics for each metric:

    * PP_sitewise (from `.state`)
    * pLDDT_CA (from `NODE_plddt_bfactor.pdb`)
    * PP_CA (from `NODE_pp_bfactor.pdb`)
    * PP−pLDDT_CA (from `NODE_ppminusplddt_bfactor.pdb`)
    * PP×pLDDT_CA (from `NODE_ppxplddt_bfactor.pdb`)

---

## 7. Notes on the metrics

* **PP (posterior probability)**

  * Derived from IQ-TREE `.state`
  * For each site, the maximum posterior probability over all states is used
  * Typically represents the confidence in the reconstructed amino acid at that site

* **pLDDT (AlphaFold)**

  * Original per-residue confidence from AlphaFold (0–100)
  * Stored in the B-factor field of the original PDB

* **PP×pLDDT**

  * Heuristic combined metric: high only if both ASR PP and pLDDT are high
  * Useful to identify residues with robust support from both evolutionary and structural perspectives

* **PP−pLDDT (scaled)**

  * Highlights disagreement between ASR confidence and AlphaFold confidence
  * After centering (−100..100), values are scaled to 0–100 solely for visualization convenience

These metrics are primarily intended for **qualitative visualization** (PyMOL, ChimeraX, etc.) and for simple region-wise averages (e.g. TM vs loop).

---

## 8. Visualization tips (optional)

Example: coloring by B-factor in PyMOL for `NODE_ppxplddt_bfactor.pdb`:

```pymol
load Node112_ppxplddt_bfactor.pdb, pp_plddt
spectrum b, blue_red, pp_plddt
```

In ChimeraX, you can use:

```chimerax
open Node112_ppxplddt_bfactor.pdb
color byattribute bfactor
```

Color schemes can be tuned to highlight high-confidence transmembrane segments versus low-confidence loops, etc.

---

## 9. Troubleshooting

* **`[ERROR] No *_summary_confidences_*.json found`**

  * Make sure you are running the script inside the AlphaFold output directory for the node (e.g., `Node112/`).
  * Check that AlphaFold completed successfully and produced `*_summary_confidences_*.json`.

* **`[ERROR] NODE not found in ALN_FASTA (exact header match)`**

  * Confirm that your alignment FASTA contains a header exactly equal to `>NodeXYZ` (no extra spaces or suffixes).
  * The script does an exact string match on the header line.

* **`[ERROR] Ungapped ASR length != structure sequence length`**

  * Indicates a mismatch between the ungapped ASR sequence and the AlphaFold structure sequence.
  * Check that:

    * the sequence used for AlphaFold exactly matches the ancestral sequence for this node
    * no extra residues (e.g. tags, signal peptides) are present in the structure that are absent in the ASR alignment.

* **Biopython / gemmi import errors**

  * Install the required Python libraries, e.g.:

    ```bash
    pip install gemmi biopython
    ```

---
# confmap (PAML + AlphaFold) pipeline

This directory contains a small pipeline to map **ancestral sequence reconstruction (ASR)** posterior probabilities (PP) from **PAML `.rst` files** onto **AlphaFold** models as B-factors, and to derive metrics such as **PP × pLDDT** and **PP − pLDDT** for structural visualization and analysis.

The workflow is designed for **PAML-based ASR on an indel-aware alignment**, combined with **AlphaFold monomer** predictions for selected ancestral nodes (e.g., Node233, Anc-SzR, Anc-HeR, etc.).

---

## 1. Files

- `run_confmap_paml.sh`  
  Main wrapper script. For a given internal node (PAML node ID, integer):
  - finds the corresponding AlphaFold model (`*_summary_confidences_*.json`, `*_model_X.cif`)
  - converts mmCIF → PDB (keeping pLDDT in the B-factor field)
  - extracts **per-alignment-column** PP values for that node from the PAML `.rst` file
  - checks alignment length consistency against a **with-gap** FASTA
  - maps PP to B-factors (per residue) using `map_confidence_to_bfactor.py`
  - generates additional PDB files with **PP − pLDDT (scaled)** and **PP × pLDDT**
  - computes CA-only statistics and writes a log file

- `extract_pp_from_paml_rst.py`  
  Python script that:
  - reads a PAML `.rst` file
  - determines the internal node order from `node #N` lines
  - extracts the PP value for the specified node **for each alignment column** from the “site / Freq / Data” block
  - optionally checks length consistency against a with-gap FASTA  
  - writes one PP value per line (alignment column)

- `map_confidence_to_bfactor.py`  
  Python script (using `gemmi` and `Bio.PDB`) that:
  - reads an AlphaFold mmCIF file
  - reads an aligned ASR sequence (with gaps) for the node of interest
  - reads per-site confidence values (PP)
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

- **PAML ASR** (`codeml` or related program), producing:
  - a PAML `.rst` file that contains a **marginal reconstruction of ancestral sequences** with per-site posterior probabilities, e.g.:

    ```text
    HeR_SzR_228_FFT_NS_1.rst
    ```

- **Indel-aware alignment**  
  An alignment including all extant and ancestral sequences, where each **PAML internal node ID** appears as a FASTA header (e.g., `>233`):

  ```text
  ASR_LGFG4_FFT_NS_1_indel_withgap.fasta
  ```

  The alignment should be the same logical alignment used to generate the PAML `.rst` file (possibly after mapping from an upstream RAxML pipeline).

- **AlphaFold** run for the ancestral sequence of a specific node (Node233, Anc-SzR, etc.), producing a directory like:

  ```text
  Node233/
    fold_YYYY_MM_DD_HH_MM_..._summary_confidences_0.json
    fold_YYYY_MM_DD_HH_MM_..._summary_confidences_1.json
    ...
    fold_YYYY_MM_DD_HH_MM_..._model_0.cif
    fold_YYYY_MM_DD_HH_MM_..._model_1.cif
    ...
    msas/
    templates/
    terms_of_use.md
  ```

  Filenames **may** contain the node label (e.g. `node233`); if so, those are preferred when choosing the best model.

---

## 3. Directory layout and basic idea

The intended usage is:

1. AlphaFold is run for an ancestral node (e.g., Node233).  
   This creates a directory such as `Node233/` containing several `*_summary_confidences_*.json` files and `*_model_X.cif` files.

2. `run_confmap_paml.sh` is executed **inside that AlphaFold output directory**, using:

   * a PAML `.rst` file
   * a with-gap alignment FASTA file
   * the PAML internal node ID (integer, e.g. `233`)

3. The script creates a subdirectory (by default `confmap/`) and writes all derived PDBs and a log file there.

Typical directory:

```text
Node233/
  fold_..._node233_summary_confidences_0.json
  fold_..._node233_model_0.cif
  ...
  msas/
  templates/
  terms_of_use.md
  confmap/
    Node233_plddt_bfactor.pdb
    Node233_pp_bfactor.pdb
    Node233_ppminusplddt_bfactor.pdb
    Node233_ppxplddt_bfactor.pdb
    Node233_aligned_indel.fasta
    Node233_pp_values.txt
    Node233_conf_CA_stats.log
```

Here, `Node233` is a label constructed from the PAML node ID (`233`) and is used consistently for output filenames and AlphaFold model selection.

---

## 4. Inputs

### 4.1 PAML ASR outputs

* `--rst <PAML_RST>`  
  PAML `.rst` file containing marginal ancestral reconstructions, e.g.:

  ```bash
  HeR_SzR_228_FFT_NS_1.rst
  ```

  The script expects a block like:

  ```text
  (1) Marginal reconstruction of ancestral sequences
  (eqn. 4 in Yang et al. 1995 Genetics 141:1641-1650).

  Prob of best state at each node, listed by site

     site   Freq   Data:
     1   ...
     2   ...
     ...
  ```

  where each site line contains a list of `X(0.xxx)` entries, one per node.

* `--node <PAML_NODE_INT>`  
  PAML internal node ID (integer), e.g.:

  ```bash
  --node 233
  ```

  This node ID must:
  - appear in the RST file as `node #233`
  - correspond to a header `>233` in the with-gap alignment FASTA

* `--withgap <WITHGAP_FASTA>`  
  Indel-aware alignment (FASTA) that includes the **PAML node ID as a header**:

  ```bash
  ASR_LGFG4_FFT_NS_1_indel_withgap.fasta
  ```

  The script expects an exact header:

  ```text
  >233
  ```

  and extracts that single sequence to `Node233_aligned_indel.fasta`.

### 4.2 AlphaFold outputs

From the AlphaFold run for the ancestral sequence of the node:

* `*_summary_confidences_*.json`  
  (optionally containing the node label, e.g. `_node233_`)

* `*_model_X.cif`  
  (corresponding AlphaFold models in mmCIF format)

The script automatically:

1. Collects all `*_summary_confidences_*.json`
2. Prefers those whose filename contains the node label (e.g. `node233`, case-insensitive)
3. Among the candidates, selects the one with the highest `ranking_score` (if present) or `plddt`
4. Extracts the model index `X` from `_summary_confidences_X.json`
5. Uses the corresponding `*_model_X.cif` as the best AlphaFold model

---

## 5. Usage

Run **inside** the AlphaFold directory for a given node:

```bash
cd Node233/

bash /path/to/run_confmap_paml.sh \
  --rst     /path/to/HeR_SzR_228_FFT_NS_1.rst \
  --node    233 \
  --withgap /path/to/ASR_LGFG4_FFT_NS_1_indel_withgap.fasta \
  --outdir  confmap_233
```

Arguments:

* `--rst`     : path to PAML `.rst` file
* `--node`    : PAML node ID (integer, e.g. `233`)
* `--withgap` : with-gap alignment FASTA containing a header `>PAML_NODE`
* `--outdir`  : output directory for confmap results (default: `confmap`)

Help:

```bash
bash run_confmap_paml.sh --help
```

---

## 6. Outputs

Within `--outdir` (default `confmap/`), the script writes:

* `Node233_plddt_bfactor.pdb`

  * PDB converted from the selected AlphaFold `.cif`
  * B-factors = original pLDDT values by AlphaFold

* `Node233_pp_bfactor.pdb`

  * PDB where B-factors = ASR posterior probability (PP) per residue
  * PP is automatically scaled in `map_confidence_to_bfactor.py`:

    * if max(PP) ≤ 1.0, treat as 0–1 and scale to 0–100
    * otherwise assume already in 0–100

* `Node233_ppminusplddt_bfactor.pdb`

  * B-factor encodes **scaled (PP − pLDDT)**:

    * raw difference is `diff = PP − pLDDT` (range approx −100 to +100)
    * then linearly scaled to `[0, 100]` via `(diff + 100)/200 × 100`

* `Node233_ppxplddt_bfactor.pdb`

  * B-factor encodes **PP × pLDDT / 100**, i.e. a simple product metric:

    * PP and pLDDT are treated as 0–100 values
    * product is `(PP * pLDDT) / 100`, constrained to 0–100

* `Node233_aligned_indel.fasta`

  * the extracted aligned sequence for this node (one entry, with gaps)
  * taken directly from `--withgap` FASTA (header `>233`)

* `Node233_pp_values.txt`

  * site-wise PP values (one per alignment column) extracted from the PAML `.rst`:

    * for each site line, the script takes the PP corresponding to the PAML node index (based on the `node #N` order)
    * no additional aggregation over states is done; values are directly the posterior probability for that node at each column

* `Node233_conf_CA_stats.log`

  * CA-only statistics for each metric:

    * `PP_sitewise`      (from `Node233_pp_values.txt`)
    * `pLDDT_CA`         (from `Node233_plddt_bfactor.pdb`)
    * `PP_CA`            (from `Node233_pp_bfactor.pdb`)
    * `PP-PLDDT_CA`      (from `Node233_ppminusplddt_bfactor.pdb`)
    * `PPxPLDDT_CA`      (from `Node233_ppxplddt_bfactor.pdb`)

---

## 7. Notes on the metrics

* **PP (posterior probability)**

  * Derived from PAML `.rst`:
    - for each site, PAML lists best states and their posterior probabilities for all nodes
    - the script extracts the probability corresponding to the chosen node ID
  * Typically represents the confidence in the ancestral residue at that site for this particular node

* **pLDDT (AlphaFold)**

  * Original per-residue confidence from AlphaFold (0–100)
  * Stored in the B-factor field of the original PDB converted from the mmCIF

* **PP×pLDDT**

  * Heuristic combined metric: high only if both PAML PP and pLDDT are high
  * Useful for identifying residues supported both by evolutionary reconstruction and by the AlphaFold model

* **PP−pLDDT (scaled)**

  * Highlights disagreement or imbalance between ASR confidence and AlphaFold confidence
  * After centering (−100..100), values are scaled to 0–100 purely for visualization convenience

These metrics are primarily intended for **qualitative visualization** (PyMOL, ChimeraX, etc.) and simple region-wise summaries (e.g., transmembrane vs extramembrane).

---

## 8. Visualization tips (optional)

Example: coloring by B-factor in PyMOL for `Node233_ppxplddt_bfactor.pdb`:

```pymol
load Node233_ppxplddt_bfactor.pdb, pp_plddt
spectrum b, blue_red, pp_plddt
```

In ChimeraX, you can use:

```chimerax
open Node233_ppxplddt_bfactor.pdb
color byattribute bfactor
```

You can apply similar commands to:

- `Node233_plddt_bfactor.pdb` (pLDDT only)
- `Node233_pp_bfactor.pdb` (PP only)
- `Node233_ppminusplddt_bfactor.pdb` (scaled PP − pLDDT)

---

## 9. Troubleshooting

* **`[ERROR] No *_summary_confidences_*.json found`**

  * Make sure you are running the script **inside** the AlphaFold output directory for the node (e.g., `Node233/`).
  * Check that AlphaFold completed successfully and produced `*_summary_confidences_*.json`.

* **`[ERROR] >PAML_NODE not found in WITHGAP_FASTA (exact header match).`**

  * Confirm that your with-gap alignment FASTA contains a header exactly equal to `>233` (no extra spaces or suffixes).
  * The script does an exact string match on the header line.

* **`[ERROR] Node N was not found in the RST file.`**

  * Check that the requested `--node` (e.g., `233`) appears as `node #233` in the `.rst` file.
  * If you changed node numbering between RAxML and PAML, ensure you are using the PAML node ID.

* **Length mismatch warnings**

  * If you see:

    ```text
    [WARN] PAML PP length (...) and with-gap alignment length (...) differ. Truncating to the shorter length.
    ```

    then:

    - confirm that the alignment used for PAML and the with-gap FASTA are consistent
    - check for trimmed regions or additional positions (e.g., masked columns) in one but not the other

* **Biopython / gemmi import errors**

  * Install the required Python libraries, e.g.:

    ```bash
    pip install gemmi biopython
    ```

---

This PAML-based pipeline is intended to be fully parallel to the IQ-TREE-based version, differing only in how the site-wise PP values are extracted (from `.rst` instead of `.state`).

---
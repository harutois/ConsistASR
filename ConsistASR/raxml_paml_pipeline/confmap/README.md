# confmap (PAML + AlphaFold) pipeline

This directory contains a small pipeline to map **ancestral sequence reconstruction (ASR)** posterior probabilities (PP) from **PAML `rst` files** onto **AlphaFold** models as B-factors, and to derive metrics such as **PP × pLDDT** and **PP − pLDDT** for structural visualization and analysis.

The workflow is designed for **PAML-based ASR and gap-containing ancestral sequences** that preserve the same alignment-column coordinates as the alignment represented in the PAML `rst` file, combined with **AlphaFold monomer** predictions for selected ancestral nodes.

The with-gap FASTA can represent either raw PAML ASR sequences or indel-aware processed sequences, provided that it preserves the same alignment-column coordinates as the PAML `rst` output and that its gap-stripped sequence matches the sequence used for AlphaFold prediction.

---

## 1. Files

- `run_confmap_paml.sh`  
  Main wrapper script. For a given internal node (PAML node ID, integer):
  - finds the corresponding AlphaFold model (`*_summary_confidences_*.json`, `*_model_X.cif`)
  - converts mmCIF → PDB (keeping pLDDT in the B-factor field)
  - extracts **per-alignment-column** PP values for that node from the PAML `rst` file
  - checks alignment length consistency against a **with-gap** FASTA
  - maps PP to B-factors (per residue) using `map_confidence_to_bfactor.py`
  - generates additional PDB files with **PP − pLDDT (scaled)** and **PP × pLDDT**
  - computes CA-only statistics and writes a log file

- `extract_pp_from_paml_rst.py`  
  Python script that:
  - reads a PAML `rst` file
  - determines the ancestral-node order from the PAML declaratio `Nodes N to M are ancestral`
  - cross-checks this order against reconstructed `node #N` sequence lines when present
  - extracts the PP value for the specified node for each alignment column from the `site / Freq / Data` block
  - verifies that each site contains the expected number of node probabilities
  - optionally checks alignment-column consistency against a with-gap FASTA
  - writes one PP value per alignment column

- `map_confidence_to_bfactor.py`  
  Python script (using `gemmi`) that:
  - reads an AlphaFold mmCIF file
  - reads an aligned ASR sequence (with gaps) for the node of interest
  - reads per-site PP values
  - maps PP and/or pLDDT to per-residue B-factors, handling gaps explicitly
  
- `extract_bfactor_from_pdb.py`  
  Auxiliary utility that extracts B-factor values from one or more PDB files. By default, it extracts CA atoms and writes a TSV table containing residue metadata and B-factor values. This is useful for exporting PP, pLDDT, PP × pLDDT, or PP − pLDDT values from the generated PDB files.

---

## 2. Dependencies

### Runtime

- Bash (tested with `bash`)
- Python `>=3.10,<3.13`
- Python libraries:
  - `gemmi` (mmCIF parsing and PDB output)
  - `biopython` (`Bio.PDB` module for simple PDB manipulation)
- One of the following mmCIF → PDB converters available on `$PATH`:
  - `phenix.cif_as_pdb`
  - `gemmi` (command-line `gemmi convert`)
  - `pdb_fromcif`

### Upstream (assumed already done)

- **PAML ASR** (`codeml` or related program), producing:
  - a PAML `rst` file that contains a **marginal reconstruction of ancestral sequences** with per-site posterior probabilities (e.g., `HeR_SzR_228_FFT_NS_1.rst`)
  - **gap-containing ancestral FASTA**
A gap-containing ancestral FASTA in which each PAML internal-node ID appears as an exact FASTA header (e.g. `>233`). This may be either an indel-aware processed ancestral FASTA or a raw PAML ASR with-gap FASTA.

  ```text
  ASR_LGFG4_FFT_NS_1_indel_withgap.fasta
  ```

The with-gap FASTA must preserve the same alignment-column coordinate system as the alignment represented in the PAML `rst` file.

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

- `extract_bfactor_from_pdb.py` uses only the Python standard library.

---

## 3. Directory layout and basic idea

The intended usage is:

1. AlphaFold is run for an ancestral node (e.g., Node233).  
   This creates a directory such as `Node233/` containing several `*_summary_confidences_*.json` files and `*_model_X.cif` files.

2. `run_confmap_paml.sh` is executed **inside that AlphaFold output directory**, using:

   * a PAML `rst` file
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
    Node233_aligned_withgap.fasta
    Node233_pp_values.txt
    Node233_conf_CA_stats.log
    Node233_bfactor_values.tsv  (optional, generated by extract_bfactor_from_pdb.py)
```

Here, `Node233` is a label constructed from the PAML node ID (`233`) and is used consistently for ConsistASR output filenames. The user-defined AlphaFold Server job name is independent of the PAML node label.

---

## 4. Inputs

### 4.1 PAML ASR outputs

* `--rst <PAML_RST>`  
  PAML `rst` file containing marginal ancestral reconstructions, e.g.:

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

 - fall within the ancestral-node range declared in the RST file (e.g. `Nodes 19 to 35 are ancestral`)
 - correspond to an exact header `>233` in the with-gap ancestral FASTA

  When reconstructed `node #N` sequence lines are present in the RST file, the helper script also verifies that their order agrees with the declared ancestral-node range.

* `--withgap <WITHGAP_FASTA>`  
  With-gap ASR alignment FASTA containing the **PAML node ID as a header**.

  This can be either:

  * an indel-aware with-gap FASTA file, e.g.  
    `<prefix>_indel_withgap.fasta`

  * a raw ASR with-gap FASTA file, e.g.  
    `<prefix>_raw_asr_withgap.fasta`

  Use the with-gap FASTA that corresponds to the gap-stripped sequence used for AlphaFold prediction.
  For example, if AlphaFold was run using `<prefix>_raw_asr_nogap.fasta`, use `<prefix>_raw_asr_withgap.fasta` here.

  The script expects an exact header `>${NODE}` (e.g., `>233`).


### 4.2 AlphaFold outputs

From the AlphaFold run for the ancestral sequence of the node:

* `*_summary_confidences_*.json`  

* `*_model_X.cif` (corresponding AlphaFold models in mmCIF format)

The script automatically:

1. Collects all `*_summary_confidences_*.json` files in the current directory.
2. Confirms that they belong to a single AlphaFold Server output prefix.
3. Checks for the corresponding `*_model_X.cif` file for each summary JSON.
4. Among valid JSON/CIF pairs, selects the model with the highest `ranking_score` (or `plddt` if `ranking_score` is unavailable).
5. Verifies downstream that the selected structure sequence exactly matches the gap-stripped ancestral sequence specified by `--withgap`.

The user-defined AlphaFold Server job name is not used to identify the ancestral node, because AlphaFold Server output filenames reflect arbitrary user-assigned job names.

---

## 5. Usage

Run **inside** the AlphaFold directory for a given node:

```bash
cd Node233/

bash /path/to/run_confmap_paml.sh \
  --rst     /path/to/HeR_SzR_228_FFT_NS_1.rst \
  --node    233 \
  --withgap /path/to/ASR_LGFG4_FFT_NS_1_indel_withgap.fasta \
  --outdir  confmap
```

Arguments:

* `--rst`     : path to PAML `rst` file
* `--node`    : PAML node ID (integer, e.g. `233`)
* `--withgap` : with-gap alignment FASTA containing a header `>PAML_NODE`
* `--outdir`  : output directory for confmap results (default: `confmap`)

Help:

```bash
bash run_confmap_paml.sh --help
```

The example above uses an indel-aware FASTA file. To map confidence values onto a structure predicted from a raw ASR sequence, replace `--withgap` with the corresponding raw with-gap FASTA file, for example:

```bash
--withgap /path/to/ASR_LGFG4_FFT_NS_1_raw_asr_withgap.fasta
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

* `Node233_aligned_withgap.fasta`

  * the extracted aligned sequence for this node (one entry, with gaps)
  * taken directly from `--withgap` FASTA (header `>233`)

* `Node233_pp_values.txt`

  * site-wise PP values (one per alignment column) extracted from the PAML `rst`:

    * the ancestral-node order is derived from the PAML declaration `Nodes N to M are ancestral`
    * for each site, the script verifies that the number of posterior probabilities matches the number of ancestral nodes
    * the PP corresponding to the requested node is then extracted
    * no additional aggregation over states is performed

* `Node233_conf_CA_stats.log`

  * summary statistics containing:

    * `PP_sitewise_0to100` — all alignment columns in the PAML site-wise output
    * `pLDDT_CA` — retained structure residues represented by CA atoms
    * `PP_CA`
    * `PP-PLDDT_CA`
    * `PPxPLDDT_CA`

  * `PP_sitewise_0to100` and the CA-based metrics may contain different numbers of positions because alignment columns corresponding to gaps are not represented as structure residues.

---

## 7. Optional: extracting B-factor values as TSV

The generated PDB files store PP, pLDDT, PP × pLDDT, or scaled PP − pLDDT values in the B-factor column.
To export these values as a TSV table, use:

```bash
python3 /path/to/extract_bfactor_from_pdb.py \
  --in Node233_pp_bfactor.pdb \
  --out Node233_pp_bfactor.tsv
```

By default, the script extracts CA atoms and writes:

`file    model   chain   resseq  icode   resname atom    bfactor`

Multiple PDB files can be processed together:

```bash
python3 /path/to/extract_bfactor_from_pdb.py \
  --in Node233_pp_bfactor.pdb Node233_plddt_bfactor.pdb Node233_ppxplddt_bfactor.pdb \
  --out Node233_bfactor_values.tsv
```

To write only B-factor values, one per line:

```bash
python3 /path/to/extract_bfactor_from_pdb.py \
  --in Node233_pp_bfactor.pdb \
  --values-only > Node233_pp_bfactor_values.txt
```

Useful options:

  * --atom CA
    Atom name to extract. Default: CA.
  * --chain A
    Optional chain ID filter.
  * --values-only
    Write only B-factor values, one per line.
  * --no-header
    Omit the TSV header line.

---

## 8. Notes on the metrics

The PDB B-factor field is used here as a convenient container for visualization; these values are not experimental crystallographic B-factors.

* **PP (posterior probability)**

  * Derived from PAML `rst`:
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
  * The raw difference ranges approximately from −100 to +100 and is linearly rescaled to 0–100 for visualization convenience.

PP×pLDDT and scaled PP−pLDDT are heuristic visualization and summary metrics. They should not be interpreted as calibrated probabilities, formal statistical confidence measures, or combinations of statistically independent quantities.

These metrics are primarily intended for **qualitative visualization** (PyMOL, ChimeraX, etc.) and simple region-wise summaries (e.g., transmembrane vs extramembrane).

---

## 9. Visualization tips (optional)

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

Color schemes can be tuned to highlight high-confidence transmembrane segments versus low-confidence loops, etc.

You can apply similar commands to:

- `Node233_plddt_bfactor.pdb` (pLDDT only)
- `Node233_pp_bfactor.pdb` (PP only)
- `Node233_ppminusplddt_bfactor.pdb` (scaled PP − pLDDT)

---

## 10. Troubleshooting

* **`[ERROR] No *_summary_confidences_*.json found`**

  * Make sure you are running the script **inside** the AlphaFold output directory for the node (e.g., `Node233/`).
  * Check that AlphaFold completed successfully and produced `*_summary_confidences_*.json`.

* **`[ERROR] Exact FASTA header >PAML_NODE not found ...`**

  * Confirm that your with-gap alignment FASTA contains a header exactly equal to `>233` (no extra spaces or suffixes).
  * The script does an exact string match on the header line.

* **`[ERROR] Node N was not found in the RST file.`**

  * Check that the requested `--node` (e.g., `233`) appears as `node #233` in the `rst` file.
  * If you changed node numbering between RAxML and PAML, ensure you are using the PAML node ID.

* **`[ERROR] PAML PP length (...) != with-gap alignment length (...)`**

  * The PAML `rst` output and the with-gap ancestral FASTA do not use the same alignment-column coordinate system.
  * Confirm that the FASTA corresponds to the same alignment used for the PAML reconstruction.
  * Do not trim or truncate either input to force the lengths to match.

* **Biopython / gemmi import errors**

  * Install the required Python libraries, e.g.:

    ```bash
    pip install gemmi biopython
    ```

* **`[ERROR] No matching ATOM records found` when using `extract_bfactor_from_pdb.py`**

  * Check that the input file is a PDB file, not mmCIF.
  * Confirm that the requested atom exists, e.g. `--atom CA`.
  * If using `--chain`, confirm that the chain ID matches the PDB file.

* **`[ERROR] Unexpected number of node probabilities at site ...`**

  * The marginal-reconstruction site row does not contain the expected
    number of posterior probabilities.
  * Confirm that the `rst` file is complete and unmodified.
  * The script stops rather than substituting a placeholder value.

* **`[ERROR] AlphaFold structure sequence does not exactly match the
  ungapped ancestral sequence`**

  * Confirm that the AlphaFold model was generated from the same ancestral sequence represented by the selected `--withgap` entry.
  * Equal sequence lengths alone are insufficient; the amino-acid sequences must match exactly.

* **`[ERROR] Multiple AlphaFold output prefixes were found`**

  * Run the wrapper in a directory containing outputs from only one AlphaFold Server job.
  * The AlphaFold job name is user-defined and is not used to infer the PAML node identity.

---

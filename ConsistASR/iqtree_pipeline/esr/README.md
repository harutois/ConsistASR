
# ESR proxy validation utilities (IQ-TREE)

This folder contains helper scripts to perform a lightweight **Extant Sequence Reconstruction (ESR)** proxy test.
The idea is to remove one extant sequence from the dataset and replace it with a synthetic proxy subtree, then
use IQ-TREE ASR on that modified dataset to see how well the internal proxy node reconstructs the held-out extant sequence.

Scripts:

- `make_esr_proxy.py` : generate ESR proxy alignment/tree (FASTA + NEWICK)
- `score_esr.py` : score ESR accuracy (global)
- `score_esr_by_region.py` : score ESR accuracy by regions

---

## Concept (what this ESR proxy does)

Given an original alignment + tree:

1. Choose a target extant tip `TARGET` (e.g. SzR4, HeR 48C12).
2. In the tree, replace `TARGET` with a proxy subtree consisting of two dummy leaves:
   - `TARGET_A` and `TARGET_B` attached with a **long branch** (`--dummy-bl`).
3. In the alignment, remove the real `TARGET` sequence and add dummy sequences for `TARGET_A` and `TARGET_B`
   (dummy residues are produced internally by `make_esr_proxy.py`).
4. Run IQ-TREE ASR on the proxy dataset.
5. Read the reconstructed residues at the **internal node labeled `TARGET`** and compare them to the original
   `TARGET` residues in the original alignment.

Key point:
- This is a **proxy validation** (not a full simulation study). It checks whether the pipeline can reconstruct
  known extant sequences from their phylogenetic context under the chosen alignment/model settings.

---

## Requirements

- Python 3.10+
- IQ-TREE v3 (`iqtree3`)
- Original MSA FASTA (gapped)
- Original tree NEWICK (with branch lengths)

---

## Step-by-step usage

### 0) Prepare inputs

Example filenames (project-specific):
- MSA: `HeR_SzR_OG.PSITM.fasta`
- Tree: `HeR_SzR_OG_QPFAMR7_PSITM.treefile`
- Target IDs:
  - SzR4: `SzR_AM_5_00977`
  - HeR 48C12: `HeR_AVZ43932`

### 1) Create ESR proxy dataset

Example (SzR4):

```bash
python make_esr_proxy.py \
  --msa HeR_SzR_OG.PSITM.fasta \
  --tree HeR_SzR_OG_QPFAMR7_PSITM.treefile \
  --target SzR_AM_5_00977 \
  --dummy-bl 500 \
  --out-prefix ESR_SzR4
```

Outputs:

* `ESR_SzR4.fasta`  (proxy alignment; contains `SzR_AM_5_00977_A` and `SzR_AM_5_00977_B`)
* `ESR_SzR4.tree`   (proxy tree; replaces the original `SzR_AM_5_00977` tip by a proxy subtree)

Notes:

Recommended branch-length handling depends on the IQ-TREE model.

* If the selected model supports fixed branch lengths, use `-blfix` and a moderate dummy branch length such as `--dummy-bl 10`.

* If branch-length fixing is not supported, for example under FreeRate models such as `+R`, do not use `-blfix`. In that case, longer dummy branches such as `--dummy-bl 100–500` may be useful, but the behavior should be checked empirically.

### 2) Run IQ-TREE ASR on the proxy dataset

Example:

```bash
iqtree3 -s ESR_SzR4.fasta \
  -m Q.pfam+R7 \
  -te ESR_SzR4.tree \
  -asr \
  -nt AUTO \
  --prefix ESR_SzR4_QPFAMR7
```

Expected output:

* `ESR_SzR4_QPFAMR7.state`  (posterior probabilities for all nodes, all sites)

Important:

* Do NOT use `-blfix` when using FreeRate models (e.g. `+R`). IQ-TREE may throw:
  `ERROR: Fixing branch lengths not supported under specified site rate model`.

### 3) Score ESR accuracy (global)

```bash
python score_esr.py \
  --orig-msa HeR_SzR_OG.PSITM.fasta \
  --target SzR_AM_5_00977 \
  --state ESR_SzR4_QPFAMR7.state
```

Output fields:

* `Comparable sites (true != gap)`: number of columns where the true target is not a gap/missing
* `Identity (%)`: fraction of comparable sites where the MAP residue equals the true residue
* `Mean posterior prob assigned to TRUE residue`: mean posterior probability assigned to the true residue

### 4) Optional: score ESR by user-defined regions

Regional scoring is optional. It is useful when meaningful sequence regions are known, such as transmembrane helices, extramembrane loops, secondary-structure elements, domains, or linkers.

Region definitions are provided as an external TSV file rather than being hard-coded in the script.

The region TSV must contain four columns:

```text
region	group	start	end
```

Coordinates are 1-based ungapped residue coordinates of the true target sequence.

Example:

```text
region	group	start	end
N-terminal EM	EM	1	2
TM1	TM	3	24
EM1	EM	25	29
TM2	TM	30	50
```

Example command:

```bash
python score_esr_by_region.py \
  --orig-msa HeR_SzR_OG.PSITM.fasta \
  --target SzR_AM_5_00977 \
  --state ESR_SzR4_QPFAMR7.state \
  --regions examples/regions/rhodopsin_SzR_regions.tsv
```

Outputs:

* GLOBAL summary
* Group-level summary, e.g. TM and EM
* Per-region summary, e.g. TM1, EM1, etc.

Example region files are provided in:

```text
examples/regions/
```

---

## Settings used in the rhodopsin example

For the microbial rhodopsin ESR proxy tests, the following settings were used:

* Model: `Q.pfam+R7`
* Branch-length fixing: not used, because `-blfix` is unsupported under `+R`
* Dummy branch length: `--dummy-bl 500` (empirically stable for this dataset/model)
* ASR run: fixed topology with `-te <ESR.tree>`

These settings are not universal recommendations. For other models or protein families, tune `--dummy-bl` and branch-length handling as needed.

---

## Troubleshooting

### 1) “Sequence TARGET_A contains only gaps or missing data”

Cause:

* If dummy sequences were encoded as missing (e.g., X) they may be treated as missing data.
In this project, make_esr_proxy.py uses fixed dummy amino acids (L and V), so this error should not occur unless the MSA length/header matching is broken.

Fix:

* Ensure dummy sequences are composed of standard amino acids (e.g. all `L` and all `V`), not `X`.
* If `make_esr_proxy.py` currently outputs `X`-only dummies, modify it to output two different real amino acids.
  (This project used `L` and `V`.)

### 2) The proxy subtree collapses / A and B become parallel in the inferred tree

This can happen because branch lengths are optimized and the proxy constraints are not fully preserved.
Mitigations:

* Increase `--dummy-bl` (e.g. 500+)
* Ensure the two dummy sequences differ (A != B) so they are not identical
* Keep topology fixed via `-te ESR_*.tree`

### 3) `SyntaxWarning: invalid escape sequence '\g'`

This warning comes from Python string escaping in a regex replacement string.
It is typically harmless. To silence it, use raw strings (prefix `r"..."`) in the relevant line.

---

## Minimal command summary

# Create proxy
```bash
python make_esr_proxy.py --msa <MSA> --tree <TREE> --target <ID> --dummy-bl 500 --out-prefix ESR_<ID>
```

# Run IQ-TREE ASR
```bash
iqtree3 -s ESR_<ID>.fasta -m Q.pfam+R7 -te ESR_<ID>.tree -asr -nt AUTO --prefix ESR_<ID>_QPFAMR7
```

# Score (global)
```bash
python score_esr.py --orig-msa <MSA> --target <ID> --state ESR_<ID>_QPFAMR7.state
```

# Score by user-defined regions
```bash
python score_esr_by_region.py \
  --orig-msa <MSA> \
  --target <ID> \
  --state ESR_<ID>_QPFAMR7.state \
  --regions examples/regions/example_regions.tsv
```

---

## Notes on interpretation

This proxy setup follows the idea of ESR: one known extant sequence is held out, represented as an internal proxy node, and reconstructed by ASR from the surrounding phylogenetic context.

* ESR proxy validation is a diagnostic baseline, not a full simulation-based estimate of deep-node ASR accuracy.
* Region-level accuracy depends on the supplied region definitions and the target sequence coordinates.
* Always report the alignment method, tree, substitution model, branch-length handling, and `--dummy-bl` value used.

---
## Background reference

The ESR strategy used in the proxy validation utilities follows the general idea described in:

- Sennett, M. A., and Theobald, D. L. Extant sequence reconstruction: the accuracy of ancestral sequence reconstructions evaluated by extant sequence cross-validation. *J. Mol. Evol.* 92, 181–206, (2024). https://doi.org/10.1007/s00239-024-10162-3

---
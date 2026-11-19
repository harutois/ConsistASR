# ConsistASR

ConsistASR is a transparent post-processing workflow for **indel-aware ancestral sequence reconstruction (ASR)** that starts from standard IQ-TREE `.state` or PAML `rst` outputs. The scripts also provides **confidence mapping** on AlphaFold predicted ancestral structures.

The repository provides:

- An **IQ-TREE–based indel-aware pipeline** (`iqtree_pipeline/indel_aware`)
- A corresponding **RAxML-NG + PAML indel-aware pipeline** (`raxml_paml_pipeline/indel_aware`)
- **Confidence mapping (“confmap”)** scripts to embed ASR posterior probability (PP) and pLDDT into B-factors of AlphaFold models
- **FASTA ID utilities** to sanitize headers before going through PHYLIP-style formats
- A simple **pre-clustering / filtering script** to reduce very large FASTA sets before tree building
- **Raw ASR FASTA extraction** from IQ-TREE `.state` and PAML `rst` outputs, with optional indel-aware processing
- **ESR proxy validation utilities** for checking whether known extant sequences can be reconstructed from phylogenetic context

The design philosophy is:

- Prefer **small, explicit scripts** over a fully “black-box” pipeline
- Each major step is a separate script with a clear interface
- Users are encouraged to check models, trees, and ASR settings at each stage

---

## Repository layout

```text
repository root
├── README.md
├── LICENSE
├── environment.yml
├── images/
└── ConsistASR/
    ├── README.md
    ├── common
    │   ├── extract_bfactor_from_pdb.py
    │   ├── msa_to_binary.py
    │   ├── raxmlng_binary_asr.py
    │   └── map_confidence_to_bfactor.py
    ├── fasta_id_tools
    │   ├── README.md
    │   ├── examples
    │   ├── fasta_rename_sequential.py
    │   ├── fasta_rename_uniprot.py
    │   ├── fasta_sanitize_id_underscore.py
    │   └── fasta_truncate_id_simple.py
    ├── precluster_tools
    │   ├── README.md
    │   └── filter_cluster_generic.sh
    ├── iqtree_pipeline
    │   ├── indel_aware
    │   │   ├── README.md
    │   │   ├── run_indel_aware_iqtree.sh
    │   │   └── state_and_indel_to_fasta.py
    │   ├── confmap
    │   │   ├── README.md
    │   │   └── run_confmap_iqtree.sh
    │   ├── examples
    │   │   ├── README.md
    │   │   ├── ASR/         # IQ-TREE ASR outputs for toy 7TM rhodopsin dataset
    │   │   ├── Tree/        # IQ-TREE tree search outputs
    │   │   ├── confmap/     # Example confmap outputs for Node10
    │   │   ├── indel_aware/ # Example indel-aware outputs for IQ-TREE pipeline
    │   │   └── toy_7tm_rhodopsin.fasta
    │   └── esr
    │       ├── README.md
    │       ├── make_esr_proxy.py
    │       ├── score_esr.py
    │       ├── score_esr_by_region.py
    │       └─ examples
    └── raxml_paml_pipeline
        ├── indel_aware
        │   ├── README.md
        │   ├── run_indel_aware_paml.sh
        │   ├── map_raxml_to_paml_nodes_from_rst.py
        │   └── paml_state_and_indel_to_fasta.py
        ├── confmap
        │   ├── README.md
        │   ├── run_confmap_paml.sh
        │   └── extract_pp_from_paml_rst.py
        └── examples
            ├── README.md
            ├── ASR/         # PAML ASR example (toy 7TM rhodopsin)
            ├── Tree/        # RAxML tree example
            ├── confmap/     # Example confmap outputs for Node28
            ├── indel_aware/ # Example indel-aware outputs for RAxML + PAML pipeline
            ├── run_paml_asr.py
            └── toy_7tm_rhodopsin.fasta
```
The IQ-TREE and PAML workflows use a common binary gap-state reconstruction framework, allowing indel-aware ancestral sequences obtained from different amino-acid ASR backends to be compared directly.

Each subdirectory (`iqtree_pipeline/indel_aware`, `iqtree_pipeline/confmap`, `iqtree_pipeline/examples`, `raxml_paml_pipeline/indel_aware`, `raxml_paml_pipeline/confmap`,`raxml_paml_pipeline/examples`, `fasta_id_tools`, `precluster_tools`) has its own `README.md` describing the interface and usage in more detail.

---

## Installation

We recommend using a dedicated conda/mamba environment.

### 1. Create and activate the environment

From the top-level `ConsistASR` directory:

```bash
conda config --set channel_priority strict
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate ConsistASR
```
The environment.yml file uses conda-forge and bioconda only
(nodefaults) and provides the recommended software versions for
reproducing the ConsistASR workflows.

This installs:

* Python (`>=3.10,<3.13`)
* IQ-TREE, RAxML-NG, PAML, MAFFT, FastTree
* seqkit, CD-HIT
* Biopython, gemmi, ETE3
* basic scientific Python stack (numpy, pandas)

Alternatively, all of the above can be installed on Conda (Bioconda + Conda Forge), e.g.:

```bash
conda create -n ConsistASR \
  -c conda-forge -c bioconda \
  "python>=3.10,<3.13" \
  iqtree raxml-ng paml mafft fasttree \
  seqkit cd-hit gemmi biopython ete3 \
  numpy pandas
```
For reproducible installation, use `environment.yml`; the manual command below is provided only as an alternative.

After activating the environment, the installed versions of the principal ConsistASR dependencies can be checked with:

```bash
conda list  | grep -E  \
'^(iqtree|raxml-ng|paml|seqkit|biopython|gemmi|ete3|numpy|pandas)[[:space:]]'
```

The principal command-line programs can also be verified with:

```bash
command -v iqtree3
command -v raxml-ng
command -v codeml
command -v seqkit
```

Each command should return the path to the executable within the active ConsistASR environment. The package versions should correspond to those specified in the provided environment.yml file.

The reference environment currently uses Python `>=3.10,<3.13`. The core ConsistASR workflows do not require Treemmer, but this upper bound also preserves compatibility with the optional Treemmer-based preclustering workflow, which relies on the `cgi` module removed from the Python standard library in Python 3.13.

### 2. Additional tools

Some external tools are **not** included in `environment.yml` and must be installed separately:

* **Treemmer**

  Used in `precluster_tools/filter_cluster_generic.sh` for tree-based reduction of large sequence sets.

  * Requires `ete3` (already included in the environment)
  * Treemmer (precluster_tool) depends on the deprecated cgi module and is incompatible with Python 3.13+. Please use Python <=3.12.
  * Install Treemmer from its GitHub repository and ensure the main script
    (e.g. `Treemmer_v0.3.py`) is reachable via the path used in the script or your `$PATH`

* **AlphaFold**

  The `confmap` scripts assume that AlphaFold has already been run for the ancestral sequence(s) of interest, and that the AlphaFold output directory contains:

  * `*_summary_confidences_*.json`
  * `*_model_X.cif`

 and may also contain `msas/`, `templates/`, and other AlphaFold output files

 The confmap workflow expects AlphaFold Server-style JSON/CIF output filenames. AlphaFold prediction itself is **not** performed by ConsistASR.

* **Visualization tools (optional)**

  * PyMOL, ChimeraX, or other molecular viewers
  * These are used to visualize B-factor–encoded PP / pLDDT / PP×pLDDT / PP − pLDDT, but are not required to run the scripts.

---

## Tested environment

The scripts have been tested on the following systems:

* **macOS** (Intel laptop)

  * macOS 13.7.8 (Ventura)
  * Conda/mamba environment: `ConsistASR` (see `environment.yml`)

* **Ubuntu (WSL2 on Windows)**

  * Ubuntu 24.04 LTS (WSL2)
  * Same `ConsistASR` conda environment

Key tool versions used during testing:

| Tool      | Version (macOS / Ubuntu)|
| --------- | ----------------------- |
| IQ-TREE   | 3.0.1 / 3.0.1           |
| RAxML-NG  | 1.2.2 / 1.2.2           |
| PAML      | 4.10.7 / 4.10.7         |
| MAFFT     | 7.525 / 7.525           |
| FastTree  | 2.2.0 / 2.2.0           |
| CD-HIT    | 4.8.1 / 4.8.1           |
| SeqKit    | 2.10.1 / 2.12.0         |
| Biopython | 1.85 / 1.86             |
| Gemmi     | 0.7.3 / 0.7.4           |
| ETE3      | 3.1.3 / 3.1.3           |
| Treemmer  | 0.3 (GitHub)            |


Other versions will probably work, but these are the versions we actually used to validate the examples in `iqtree_pipeline/examples` and `raxml_paml_pipeline/examples`.

---

In ConsistASR, “raw ASR” refers to ancestral amino-acid states directly extracted from IQ-TREE or PAML outputs before binary gap-state reconstruction and indel-aware processing.

## Overview of the pipelines

### 1. IQ-TREE indel-aware pipeline

`iqtree_pipeline/indel_aware/` provides a small wrapper to extract raw ancestral sequences from IQ-TREE `.state` files and to generate **indel-aware ancestral sequences** using an independent binary gap-state reconstruction.

* Input:

  * Amino-acid MSA (`--msa`)
  * IQ-TREE topology (`--tree`)
  * IQ-TREE ASR `.state` file (`--state`)
  * Outgroup taxa list (`--outgroup`)

* Main script: `run_indel_aware_iqtree.sh`

  * Converts the MSA to a binary residue-presence/residue-absence matrix using `common/msa_to_binary.py`
  * Uses RAxML-NG `--evaluate` with `BIN+G` to re-estimate branch lengths and binary-model parameters while retaining the amino-acid tree topology
  * Performs binary ancestral-state reconstruction with RAxML-NG `--ancestral`
  * Validates the RAxML-NG ancestral node labels and alignment-site coordinates against the IQ-TREE `.state` file
  * Integrates the amino-acid and binary gap-state reconstructions to generate indel-aware ancestral sequences

* Principal outputs:

  * raw ancestral FASTA files, with and without alignment gaps
  * indel-aware ancestral FASTA files, with and without reconstructed ancestral gaps
  * `*_binary_posterior.tsv`, containing posterior probabilities of residue absence and presence
  * `*_site_posterior.tsv`, containing site-wise integrated amino-acid and binary-state information
  * `*_indel_nogap_posterior.tsv`, containing posterior information corresponding directly to residues retained in the final gap-stripped ancestral sequences

See [`iqtree_pipeline/indel_aware/README.md`](iqtree_pipeline/indel_aware/README.md) and the toy example in `iqtree_pipeline/examples`.

---

### 2. IQ-TREE confmap (PP / pLDDT mapping)

`iqtree_pipeline/confmap/` provides scripts to **map ASR PP and AlphaFold pLDDT to B-factors**:

* `run_confmap_iqtree.sh`

  * Selects the best AlphaFold model from `*_summary_confidences_*.json`
  * Converts mmCIF to PDB (keeping pLDDT in B-factors)
  * Extracts site-wise PP from IQ-TREE `.state`
  * Maps PP to B-factors using `map_confidence_to_bfactor.py`
  
* Generates:

    * PDB with pLDDT
    * PDB with PP
    * PDB with scaled PP − pLDDT
    * PDB with PP × pLDDT
    * site-wise PP statistics and CA-based structure statistics

These outputs can be visualized in PyMOL or ChimeraX.
See [`iqtree_pipeline/confmap/README.md`](iqtree_pipeline/confmap/README.md)

The shared helper `common/extract_bfactor_from_pdb.py` can be used to extract residue-wise B-factor values from the generated PDB files for downstream summary or plotting.

---

### 3. RAxML-NG + PAML pipeline

`raxml_paml_pipeline/indel_aware/` and `raxml_paml_pipeline/confmap/` implement an analogous workflow for **PAML-based amino-acid ASR**. The indel-aware pipeline extracts raw ancestral sequences from PAML `rst` files and applies the same RAxML-NG-based binary gap-state reconstruction used by the IQ-TREE workflow.

* Main script: `run_indel_aware_paml.sh`

  * Uses a PAML `rst` file and the corresponding RAxML-NG phylogeny
  * Extracts raw ancestral FASTA directly from the PAML `rst` file
  * Converts the MSA to a binary residue-presence/residue-absence matrix
  * Uses RAxML-NG `--evaluate` and `--ancestral` for binary gap-state reconstruction
  * Extracts the PAML internal-node tree from `rst`
  * Maps the non-root PAML ancestral nodes to the corresponding RAxML-NG internal vertices
  * Excludes the explicit PAML root from the final indel-aware FASTA files while retaining it in the raw PAML ASR output

* Principal outputs:

  * raw ancestral FASTA files, with and without alignment gaps
  * indel-aware ancestral FASTA files, with and without reconstructed ancestral gaps
  * `*_binary_posterior.tsv`
  * `*_site_posterior.tsv`
  * `*_indel_nogap_posterior.tsv`
  * `*_node_map.tsv`, recording non-root PAML-to-RAxML-NG node correspondence

* `raxml_paml_pipeline/confmap`:

  * `run_confmap_paml.sh` reads a PAML `rst` and a with-gap ancestral FASTA
  * Extracts posterior probability per alignment column for a selected PAML node
  * Maps posterior probability and AlphaFold pLDDT to B-factors
  * Outputs PDB files and summary statistics for PP, pLDDT, PP − pLDDT, and PP × pLDDT

See [`raxml_paml_pipeline/indel_aware/README.md`](raxml_paml_pipeline/indel_aware/README.md) and [`raxml_paml_pipeline/confmap/README.md`](raxml_paml_pipeline/confmap/README.md).

---

### 4. FASTA ID utilities

`fasta_id_tools/` contains small helper scripts to avoid PHYLIP/RAxML troubles with spaces, pipes, and special characters in FASTA headers, and to keep mapping tables from original IDs to renamed IDs.

* `fasta_rename_sequential.py`

  * Renames all sequences to `seq0001`, `seq0002`, …
  * Writes a two-column mapping table (old ID → new ID)

* `fasta_truncate_id_simple.py`

  * Truncates each ID at the first non-alphanumeric/underscore character
  * E.g.

    * `AM_5S_00009|AM_5S-c1__...` → `AM_5S_00009`
    * `QBQ84358.1 schizorhodopsin [...]` → `QBQ84358`

* `fasta_sanitize_id_underscore.py`

  * Replaces spaces and problematic symbols with underscores
  * E.g.

    * `SAMEA 2622822_312577` → `SAMEA_2622822_312577`

* `fasta_rename_uniprot.py`

  * Renames UniProt-style FASTA headers (`sp|...|...` or `tr|...|...`) to shorter IDs
  * Supports ID styles based on accession, entry name, or sequential numbering
  * Writes an extended mapping table containing UniProt metadata such as accession, entry name, organism, taxonomy ID, gene name, PE, and SV
  * Supports an optional prefix for renamed IDs

Each script has its own usage examples in [`fasta_id_tools/README.md`](fasta_id_tools/README.md), including a small UniProt header-renaming example.

---

### 5. Pre-cluster tools

`precluster_tools/filter_cluster_generic.sh` implements a **generic pre-filtering & clustering workflow**:

* Cleans sequences (remove non-amino-acid characters)
* Applies a configurable minimum-length filter (default: ≥ 50 aa)
* Filters by presence of one or more residues (e.g., K; configurable)
* Removes exact duplicates
* Runs CD-HIT at several identity thresholds
* Optionally aligns and reduces the dataset with MAFFT + FastTree + Treemmer
* Writes a small report summarizing the reduction

This is intended as a **pre-processing step** to obtain a reasonable number of representative sequences before serious ASR.

See [`precluster_tools/README.md`](precluster_tools/README.md) for details.

### 6. ESR proxy validation utilities

`iqtree_pipeline/esr/` provides helper scripts for extant sequence reconstruction (ESR) proxy validation using IQ-TREE.

See [`iqtree_pipeline/esr/README.md`](iqtree_pipeline/esr/README.md) for details.

---

## Toy example: 7TM rhodopsin dataset

Both the IQ-TREE and RAxML + PAML workflows include a small 7TM rhodopsin toy dataset with bundled intermediate and output files.

- **IQ-TREE workflow:** [`iqtree_pipeline/examples/README.md`](iqtree_pipeline/examples/README.md)
- **RAxML + PAML workflow:** [`raxml_paml_pipeline/examples/README.md`](raxml_paml_pipeline/examples/README.md)

The example READMEs provide step-by-step commands for tree inference, amino-acid ASR, indel-aware processing, and confidence mapping, together with descriptions of the expected outputs.

The toy dataset is intended to demonstrate and reproduce the ConsistASR workflow rather than to serve as a biological benchmark or as a universally recommended phylogenetic analysis protocol.

---

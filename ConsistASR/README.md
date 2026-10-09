# ConsistASR

ConsistASR is a transparent post-processing workflow for **indel-aware ancestral sequence reconstruction (ASR)** that starts from standard IQ-TREE `.state` or PAML `rst` outputs. The scripts also provides **confidence mapping** on AlphaFold predicted ancestral structures.

The repository provides:

- An **IQ-TREE–based indel-aware pipeline** (`iqtree_pipeline/indel_aware`)
- A corresponding **RAxML + PAML indel-aware pipeline** (`raxml_paml_pipeline/indel_aware`)
- **Confidence mapping (“confmap”)** scripts to embed ASR posterior probability (PP) and pLDDT into B-factors of AlphaFold models
- **FASTA ID utilities** to sanitize headers before going through PHYLIP-style formats
- A simple **pre-clustering / filtering script** to reduce very large FASTA sets before tree building
- **Raw ASR FASTA extraction** IQ-TREE `.state` and PAML `rst` outputs, with optional indel-aware processing
- **ESR proxy validation utilities** for checking whether known extant sequences can be reconstructed from phylogenetic context

The design philosophy is:

- Prefer **small, explicit scripts** over a fully “black-box” pipeline
- Each major step is a separate script with a clear interface
- Users are encouraged to check models, trees, and ASR settings at each stage

---

## Repository layout

```text
ConsistASR
├── README.md
├── environment.yml
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
│   │   ├── map_raxml_to_iqtree_nodes.py
│   │   ├── msa_to_binary.py
│   │   └── state_and_indel_to_fasta.py
│   ├── confmap
│   │   ├── README.md
│   │   ├── run_confmap_iqtree.sh
│   │   ├── extract_bfactor_from_pdb.py
│   │   └── map_confidence_to_bfactor.py
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
    │   ├── msa_to_binary.py
    │   └── paml_state_and_indel_to_fasta.py
    ├── confmap
    │   ├── README.md
    │   ├── run_confmap_paml.sh
    │   ├── extract_bfactor_from_pdb.py
    │   ├── extract_pp_from_paml_rst.py
    │   └─ map_confidence_to_bfactor.py
    └── examples
        ├── README.md
        ├── ASR/         # PAML ASR example (toy 7TM rhodopsin)
        ├── Tree/        # RAxML tree example
        ├── confmap/     # Example confmap outputs for Node28
        ├── indel_aware/ # Example indel-aware outputs for RAxML + PAML pipeline
        ├── run_paml_asr.py
        └── toy_7tm_rhodopsin.fasta

```

Each subdirectory (`iqtree_pipeline/indel_aware`, `iqtree_pipeline/confmap`, `iqtree_pipeline/examples`, `raxml_paml_pipeline/indel_aware`, `raxml_paml_pipeline/confmap`,`raxml_paml_pipeline/examples`, `fasta_id_tools`, `precluster_tools`) has its own `README.md` describing the interface and usage in more detail.

---

## Installation

We recommend using a dedicated conda/mamba environment.

### 1. Create and activate the environment

From the top-level `ConsistASR` directory:

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate ConsistASR
```

This installs:

* Python (`>=3.10,<3.13`)
* IQ-TREE, RAxML-HPC, RAxML-NG, PAML, MAFFT, FastTree
* seqkit, CD-HIT
* Biopython, gemmi, ETE3
* basic scientific Python stack (numpy, pandas)

Alternatively, all of the above can be installed on Conda (Bioconda + Conda Forge), e.g.:

```bash
conda create -n ConsistASR -c conda-forge -c bioconda \
  "python>=3.10,<3.13" \
  iqtree raxml raxml-ng paml mafft fasttree seqkit cd-hit gemmi biopython ete3
```

The core Python scripts are not expected to depend strongly on the Python minor version, but the environment pins Python to `<3.13` because Treemmer depends on the deprecated `cgi` module, which was removed in Python 3.13.

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
| RAxML-HPC | 8.2.13* / 8.2.13         |
| PAML      | 4.10.7 / 4.10.7         |
| MAFFT     | 7.525 / 7.525           |
| FastTree  | 2.2.0 / 2.2.0           |
| CD-HIT    | 4.8.1 / 4.8.1           |
| SeqKit    | 2.10.1 / 2.12.0         |
| Biopython | 1.85 / 1.86             |
| Gemmi     | 0.7.3 / 0.7.4           |
| ETE3      | 3.1.3 / 3.1.3           |
| Treemmer  | 0.3 (GitHub)            |

\* The Bioconda raxml 8.2.13 package reports "RAxML version 8.2.12" in the executable banner on the tested macOS environment.

Other versions will probably work, but these are the versions we actually used to validate the examples in `iqtree_pipeline/examples` and `raxml_paml_pipeline/examples`.

---

In ConsistASR, “raw ASR” refers to ancestral amino-acid states directly extracted from IQ-TREE or PAML outputs before binary gap-state reconstruction and indel-aware processing.

## Overview of the pipelines

### 1. IQ-TREE indel-aware pipeline

`iqtree_pipeline/indel_aware/` provides a small wrapper to extract raw ancestral sequences from IQ-TREE `.state` files and to generate **indel-aware ASR** sequences using binary gap-state reconstruction.

* Input:

  * Amino-acid MSA (`--msa`)
  * IQ-TREE topology (`--tree`)
  * IQ-TREE ASR `.state` file (`--state`)
  * Outgroup taxa list (`--outgroup`)

* Main script: `run_indel_aware_iqtree.sh`

  * Converts MSA to a binary (0/1) gap matrix (`msa_to_binary.py`)
  * Evaluates the topology with RAxML-NG (BIN+G)
  * Runs RAxML-HPC binary gap-state reconstruction
  * Maps RAxML internal node IDs back to IQ-TREE node labels (`map_raxml_to_iqtree_nodes.py`)
  
* Produces:

    * raw ancestral FASTA extracted directly from the IQ-TREE `.state` file
    * indel-aware ancestral FASTA retaining alignment gaps
    * gap-stripped indel-aware ancestral FASTA for AlphaFold and downstream analyses
    * a node mapping table

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

The directory also includes `extract_bfactor_from_pdb.py`, which can extract residue-wise B-factor values from the generated PDB files for downstream summary or plotting.

---

### 3. RAxML + PAML pipeline

`raxml_paml_pipeline/indel_aware/` and `raxml_paml_pipeline/confmap/` implement an analogous workflow for **PAML-based ASR**. The indel-aware pipeline can extract raw ancestral sequences from PAML `rst` files and generate indel-aware ASR sequences using binary gap-state reconstruction.

* Main script: `run_indel_aware_paml.sh`:

  * Uses a PAML `rst` file and the RAxML tree used for the PAML run
  * Extracts raw ancestral FASTA directly from the PAML `rst` file
  * Performs binary gap-state reconstruction with RAxML-HPC
  * Extracts the PAML internal-node tree from `rst`
  * Maps RAxML internal nodes to PAML node IDs (`map_raxml_to_paml_nodes_from_rst.py`)
  
* Produces:

    * raw ancestral FASTA extracted directly from the PAML `rst` file
    * indel-aware ancestral FASTA retaining alignment gaps
    * gap-stripped indel-aware ancestral FASTA for AlphaFold and downstream analyses
    * a node mapping table

* `raxml_paml_pipeline/confmap`:

  * `run_confmap_paml.sh` reads a PAML `rst` and a with-gap ancestral FASTA
  * Extracts PP per alignment column for a given PAML node ID
  * Maps PP and AlphaFold pLDDT to B-factors (as in the IQ-TREE version)
  * Outputs PDBs for PP, PP − pLDDT, PP × pLDDT and a statistics log

See [`raxml_paml_pipeline/indel_aware/README.md`](raxml_paml_pipeline/indel_aware/README.md), [`raxml_paml_pipeline/confmap/README.md`](raxml_paml_pipeline/confmap/README.md) and `raxml_paml_pipeline/examples`.

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

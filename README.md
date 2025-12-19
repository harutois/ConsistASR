# ConsistASR

ConsistASR is a collection of small, modular scripts for **indel-aware ancestral sequence reconstruction (ASR)** and **confidence mapping** on AlphaFold models.

The main focus is on 7-transmembrane (7TM) microbial rhodopsins, but the workflow is in principle applicable to other protein families.

The repository provides:

- An **IQ-TREE–based indel-aware pipeline ("indel\_aware\_iqtree")**
- A corresponding **RAxML + PAML indel-aware pipeline ("indel\_aware\_paml")**

![indel_aware](images/indel_aware.png)

- **Confidence mapping (“confmap”)** scripts to embed ASR posterior probability (PP), PPxpLDDT and PP-pLDDT into B-factors of AlphaFold models

![confmap](images/confmap.png)

- **FASTA ID utilities** to sanitize headers before going through PHYLIP-style formats
- A simple **pre-clustering / filtering script** to reduce very large FASTA sets before tree building

---

## Installation (quick)

From the top-level `ConsistASR` directory:

```bash
git clone https://github.com/yourname/ConsistASR.git
cd ConsistASR
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate ConsistASR
```

For detailed installation, pipeline descriptions, and examples, please see the main documentation:

👉 [`ConsistASR/README.md`](ConsistASR/README.md)

---

## Citation

If you use ConsistASR in your work, please cite:

* **[Manuscript title, authors, journal, year]** – *TBD*

---

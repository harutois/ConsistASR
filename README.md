# ConsistASR

ConsistASR is a collection of small, modular scripts for **indel-aware ancestral sequence reconstruction (ASR)** and **confidence mapping** on AlphaFold models.

The main focus is on 7-transmembrane (7TM) microbial rhodopsins, but the workflow is in principle applicable to other protein families.

The repository provides:

- An **IQ-TREE–based indel-aware pipeline ("indel\_aware\_iqtree")**
- A corresponding **RAxML + PAML indel-aware pipeline ("indel\_aware\_paml")**

![indel_aware](images/indel_aware.png)

- **Confidence mapping (“confmap”)** scripts to embed ASR posterior probability (PP), PP x pLDDT and PP - pLDDT into B-factors of AlphaFold models

![confmap](images/confmap.png)

- **FASTA ID utilities** to sanitize headers before going through PHYLIP-style formats
- A simple **pre-clustering / filtering script** to reduce very large FASTA sets before tree building

---

## Installation (quick)

From the top-level `ConsistASR` directory:

```bash
git clone https://github.com/harutois/ConsistASR.git
cd ConsistASR
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate ConsistASR
```

For detailed installation, pipeline descriptions, and examples, please see the main documentation:

👉 [`ConsistASR/README.md`](ConsistASR/README.md)

---

## Reference

- Aadland, K., Pugh, C., Kolaczkowski, B. (2019). High-Throughput Reconstruction of Ancestral Protein Sequence, Structure, and Molecular Function. In: Sikosek, T. (eds) Computational Methods in Protein Evolution. Methods in Molecular Biology, vol 1851. Humana Press, New York, NY. https://doi.org/10.1007/978-1-4939-8736-8_8

---

## Citation

If you use ConsistASR in your work, please cite:

* Ishikawa, H. and Mizutani, Y. A Structure-guided, Indel-aware Framework for Ancestral Reconstruction of Full-length Seven-transmembrane Proteins. *submitted*

---

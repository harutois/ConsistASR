# Generic filter / cluster / Treemmer pipeline

This directory contains a small **generic preprocessing pipeline** for large amino-acid FASTA collections.  
The script `filter_cluster_generic.sh` takes one (or optionally two) FASTA files, applies several filters, clusters the sequences with CD-HIT, reduces the dataset using Treemmer, and finally builds a simple FastTree phylogeny.

It is intended as a quick **screening / reduction workflow** before more computationally expensive analyses.

---

## 1. Files

- `filter_cluster_generic.sh`  
  Main wrapper script. For one or two input FASTA files, it:

  1. merges the inputs (if two FASTA files are given)  
  2. cleans sequences (removes any character that is not `A–Z`/`a–z`)  
  3. applies a minimum length filter (default: `>= 200` aa)  
  4. optionally keeps only sequences that contain at least one of the specified residues (e.g. `K`)  
  5. removes exact duplicate sequences (sequence-only, ignores IDs)  
  6. runs CD-HIT at 99%, 97%, 95%, and 90% identity  
  7. aligns the 95% clustered set with MAFFT  
  8. runs Treemmer on the FastTree built from the 95% set, trimming down to `N` sequences (default: `N=300`)  
  9. builds a final FastTree on the reduced alignment  
  10. writes a simple text report

- `report.txt` (generated in the output directory)  
  Summary of key counts and file paths.

---

## 2. Dependencies

Runtime:

- Bash (tested with `bash`)
- [`seqkit`](https://bioinf.shenwei.me/seqkit/)
- [`cd-hit`](http://weizhong-lab.ucsd.edu/cd-hit/)
- [`mafft`](https://mafft.cbrc.jp/alignment/software/)
- [`FastTree`](http://www.microbesonline.org/fasttree/)
- Python 3
- [`Treemmer`](https://github.com/fmenardo/Treemmer) (e.g. `Treemmer_v0.3.py`)

The path to the Treemmer script is controlled by:

- environment variable `TREEMMER_PY`, or  
- the `--treemmer` option

---

## 3. Basic usage

### 3.1 Single FASTA as input

```bash
bash filter_cluster_generic.sh \
  --in input.fasta \
  --outdir cluster_work
```

This will:

* use `input.fasta` as the only input
* output all intermediate and final files to `cluster_work/`
* apply default settings:

  * minimum length: `>= 200` aa
  * no residue filter (all amino acids allowed)
  * Treemmer target: `N = 300` sequences

### 3.2 Merge two FASTA files

```bash
bash filter_cluster_generic.sh \
  --in   mgnify.fasta \
  --extra known.fasta \
  --outdir work_merged
```

`mgnify.fasta` and `known.fasta` are concatenated and then processed together.

---

## 4. Parameters

* `--in PATH`
  **Required.** Main input FASTA file (amino acids).

* `--extra PATH`
  Optional second FASTA file. If provided, it is concatenated with `--in`.

* `--outdir DIR`
  Output directory (default: `cluster_work`). The directory will be created if it does not exist.

* `--minlen INT`
  Minimum sequence length (aa) for the length filter.
  Default: `200`.

* `--aa-filter STR`
  Optional residue filter.
  If provided, only sequences that contain at least **one** of the characters in `STR` are kept.
  Examples:

  * `--aa-filter K`  → keep sequences with at least one Lys
  * `--aa-filter KR` → keep sequences with at least one Lys **or** Arg

  If omitted or an empty string, this step is **skipped**.

* `--target INT`
  Target number of sequences for Treemmer.
  Default: `300`.

* `--threads INT`
  Number of threads for CD-HIT and FastTree.
  Default: `8`.

* `--treemmer PATH`
  Path to the Treemmer Python script.
  If not given, the script uses:

  ```bash
  TREEMMER_PY="\$HOME/software/Treemmer/Treemmer_v0.3.py"
  ```

---

## 5. Outputs

Within `--outdir` (default `cluster_work/`), the following files are produced:

* **Filtered FASTA and clustering**

  * `all_raw.fasta`
    Simple concatenation of the input FASTA(s).

  * `all_cleaned.fasta`
    Same sequences but with any character other than `A–Z`/`a–z` removed from sequence lines.

  * `seq_lenXXX.fasta`
    Sequences after length filter (XXX = `--minlen`).

  * `seq_filtered.fasta`
    Final filtered set (after optional residue filter and duplicate removal).
    This is the input to the CD-HIT pipeline.

  * `seq_c99.fasta`, `seq_c97.fasta`, `seq_c95.fasta`, `seq_c90.fasta`
    Representative sequences after CD-HIT clustering at 99%, 97%, 95% and 90% identity.

* **Alignment and Treemmer**

  * `seq_c95.aln.fasta`
    MAFFT alignment of `seq_c95.fasta`.

  * `seq_c95.nwk`
    FastTree tree inferred from `seq_c95.aln.fasta` (WAG+Gamma).

  * `seq_c95.nwk_trimmed_list_X_N`
    Treemmer kept-list file (for target `N`; e.g., `N=300`).
    (Filename assumes Treemmer v0.3 style.)

  * `seq_N_reduced.aln.fasta`
    Alignment of the Treemmer-reduced set (or the full 95% set if Treemmer output is not found).

* **Final tree and report**

  * `seq_N.tree`
    Final FastTree (WAG+Gamma) inferred from `seq_N_reduced.aln.fasta`.

  * `seq_N_ids.txt`
    List of sequence IDs present in the final alignment / tree.

  * `report.txt`
    Summary of:

    * input files
    * filters used (`minlen`, `aa-filter`)
    * sequence counts after each CD-HIT step
    * Treemmer target `N`
    * path to final IDs and final tree

---

## 6. Typical use cases

* **Initial screening of large environmental datasets**

  Quickly reduce tens of thousands of hits to a few hundred representative sequences for:

  * manual inspection
  * more expensive ML / Bayesian phylogenetics
  * structure prediction or ASR pipelines

* **Combining “known” reference sequences with new database hits**

  Use `--extra` to merge curated reference sequences with new hits before filtering and clustering.

---

## 7. Notes and caveats

* All non-letter characters in sequence lines are removed (`[^a-zA-Z]`), so:

  * gaps (`-`), stop codons (`*`), and numbers are stripped
  * this is appropriate for amino-acid sequences but not for codon-level data

* Residue filter (`--aa-filter`) uses a simple “contains any of these characters” logic and does not consider position or motif context.

* Treemmer kept-list file naming is assumed to be compatible with **Treemmer v0.3**.
  If your version uses different names, you may need to adjust the `KEPT_LIST` pattern in the script.

---

## 8. Troubleshooting

* **`[ERROR] Treemmer script not found`**
  Check `TREEMMER_PY` or pass `--treemmer /path/to/Treemmer_v0.3.py`.

* **Very few sequences after filters**

  * Decrease `--minlen`
  * Remove `--aa-filter` or relax it (e.g. `KR` instead of `K` only)

* **CD-HIT / FastTree not found**
  Make sure the executables are available on your `\$PATH`.

---
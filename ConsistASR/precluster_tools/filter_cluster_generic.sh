#!/usr/bin/env bash
# ============================================================
# filter_cluster_generic.sh
#
# Generic filter → cluster → alignment → Treemmer → tree workflow
# for amino-acid FASTA files.
#
# Steps:
#   1) Merge one or two input FASTA files
#   2) Clean sequences (remove gaps / non-letter characters)
#   3) Length filter (default: >= 200 aa)
#   4) Optional residue filter (e.g. keep sequences containing K)
#   5) Remove exact duplicates (sequence-only)
#   6) CD-HIT clustering (99% → 97% → 95% → 90%)
#   7) MAFFT alignment of 95% clustered set
#   8) Treemmer down to N sequences (default: 300)
#   9) Final FastTree (WAG+Gamma) and simple report
#
# Dependencies:
#   - seqkit
#   - cd-hit
#   - mafft
#   - FastTree
#   - Python 3
#   - Treemmer (Python script)
#
# Usage examples:
#
#   # Single input file, default settings (minlen=200, no AA filter, target=300)
#   bash filter_cluster_generic.sh \
#       --in input.fasta \
#       --outdir cluster_work
#
#   # Merge two FASTA files, apply K filter, custom min length and target N
#   bash filter_cluster_generic.sh \
#       --in  mgnify.fasta \
#       --extra known.fasta \
#       --minlen 180 \
#       --aa-filter K \
#       --target 300 \
#       --outdir work_helios_like
#
# ============================================================

set -euo pipefail

# -------------------------
# Default parameters
# -------------------------
INPUT1=""
INPUT2=""
OUTDIR="cluster_work"
MINLEN=200          # minimum sequence length (aa)
AA_FILTER=""        # e.g. "K" or "KR"; if empty, residue filter is skipped
TARGET_N=300        # target N for Treemmer
THREADS=8           # threads for CD-HIT / FastTree
TREEMMER_PY="${TREEMMER_PY:-$HOME/software/Treemmer/Treemmer_v0.3.py}"

usage() {
    cat <<EOF
Usage:
  $0 --in INPUT.fasta [options]

Required:
  --in    PATH   Input FASTA file (amino acids)

Optional:
  --extra PATH   Second FASTA file to merge with --in
  --outdir DIR   Output directory (default: cluster_work)
  --minlen INT   Min length (aa) for length filter (default: 200)
  --aa-filter STR
                 Keep only sequences that contain at least one of these residues.
                 Example: "K" or "KR". If omitted or empty, this step is skipped.
  --target INT   Target number of sequences for Treemmer (default: 300)
  --threads INT  Number of threads for CD-HIT / FastTree (default: 8)
  --treemmer PATH
                 Path to Treemmer Python script (default: \$TREEMMER_PY)
  -h, --help     Show this help and exit

Examples:
  $0 --in input.fasta --minlen 200 --aa-filter K --target 300
EOF
}

# -------------------------
# Parse arguments
# -------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            INPUT1="$2"; shift 2;;
        --extra)
            INPUT2="$2"; shift 2;;
        --outdir)
            OUTDIR="$2"; shift 2;;
        --minlen)
            MINLEN="$2"; shift 2;;
        --aa-filter)
            AA_FILTER="$2"; shift 2;;
        --target)
            TARGET_N="$2"; shift 2;;
        --threads)
            THREADS="$2"; shift 2;;
        --treemmer)
            TREEMMER_PY="$2"; shift 2;;
        -h|--help)
            usage; exit 0;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage
            exit 1;;
    esac
done

# -------------------------
# Basic checks
# -------------------------
if [[ -z "${INPUT1}" ]]; then
    echo "[ERROR] --in is required." >&2
    usage
    exit 1
fi
if [[ ! -f "${INPUT1}" ]]; then
    echo "[ERROR] Input file not found: ${INPUT1}" >&2
    exit 1
fi
if [[ -n "${INPUT2}" ]] && [[ ! -f "${INPUT2}" ]]; then
    echo "[ERROR] Extra file not found: ${INPUT2}" >&2
    exit 1
fi
if [[ ! -f "${TREEMMER_PY}" ]]; then
    echo "[ERROR] Treemmer script not found: ${TREEMMER_PY}" >&2
    echo "        Set TREEMMER_PY env or use --treemmer PATH" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

echo "[INFO] Input1:  ${INPUT1}"
if [[ -n "${INPUT2}" ]]; then
    echo "[INFO] Input2:  ${INPUT2}"
fi
echo "[INFO] Outdir:  ${OUTDIR}"
echo "[INFO] Minlen:  ${MINLEN} aa"
if [[ -n "${AA_FILTER}" ]]; then
    echo "[INFO] AA filter: at least one of [${AA_FILTER}]"
else
    echo "[INFO] AA filter: (none; residue filter will be skipped)"
fi
echo "[INFO] Treemmer target N: ${TARGET_N}"
echo "[INFO] Treemmer script:   ${TREEMMER_PY}"
echo "[INFO] Threads:           ${THREADS}"

# -------------------------
# 1. Merge FASTA(s) → raw.fasta
# -------------------------
RAW="${OUTDIR}/all_raw.fasta"
if [[ -n "${INPUT2}" ]]; then
    cat "${INPUT1}" "${INPUT2}" > "${RAW}"
else
    cp "${INPUT1}" "${RAW}"
fi

# -------------------------
# 1b. Clean sequences: remove non-letter chars
# -------------------------
echo "[INFO] Cleaning sequences (remove non A–Z / a–z)..."
CLEANED="${OUTDIR}/all_cleaned.fasta"
awk '/^>/ { print; next }
     !/^>/ {
         gsub(/[^a-zA-Z]/, "", $0);
         print
     }' "${RAW}" > "${CLEANED}"

# -------------------------
# 2. Length filter (>= MINLEN aa)
# -------------------------
echo "[INFO] Length filter: >= ${MINLEN} aa"
IDS_LEN="${OUTDIR}/ids_len${MINLEN}.txt"

seqkit fx2tab --name --length "${CLEANED}" \
  | awk -v minlen="${MINLEN}" '$2 >= minlen { print $1 }' \
  > "${IDS_LEN}"

LEN_FASTA="${OUTDIR}/seq_len${MINLEN}.fasta"
seqkit grep -f "${IDS_LEN}" "${CLEANED}" > "${LEN_FASTA}"

CURRENT="${LEN_FASTA}"

# -------------------------
# 3. Optional residue filter
# -------------------------
if [[ -n "${AA_FILTER}" ]]; then
    echo "[INFO] Residue filter: keep sequences containing at least one of [${AA_FILTER}]"

    # e.g. "KR" -> character class "[KR]"
    AA_CLASS=$(echo "${AA_FILTER}" | tr -d ' ,')
    PATTERN="[${AA_CLASS}]"

    AA_FASTA="${OUTDIR}/seq_len${MINLEN}_aa${AA_CLASS}.fasta"
    seqkit grep -s -r -p "${PATTERN}" "${CURRENT}" > "${AA_FASTA}"

    CURRENT="${AA_FASTA}"
else
    echo "[INFO] Residue filter skipped."
fi

# -------------------------
# 4. Remove exact duplicates (sequence-only)
# -------------------------
echo "[INFO] Removing exact duplicates (sequence-only)..."
RMDUP_FASTA="${OUTDIR}/seq_filtered_rmdup.fasta"
seqkit rmdup -s -i "${CURRENT}" > "${RMDUP_FASTA}"

# -------------------------
# 5. Final filtered input
# -------------------------
FILTERED_FASTA="${OUTDIR}/seq_filtered.fasta"
cp "${RMDUP_FASTA}" "${FILTERED_FASTA}"

echo "[INFO] Filtered sequences written to: ${FILTERED_FASTA}"
echo "[INFO] Number of sequences after filters: $(grep -c '^>' "${FILTERED_FASTA}")"

# -------------------------
# 6. CD-HIT clustering
# -------------------------
echo "[INFO] Running CD-HIT at 99%, 97%, 95%, 90%..."

C99="${OUTDIR}/seq_c99.fasta"
C97="${OUTDIR}/seq_c97.fasta"
C95="${OUTDIR}/seq_c95.fasta"
C90="${OUTDIR}/seq_c90.fasta"

cd-hit -i "${FILTERED_FASTA}" -o "${C99}" -c 0.99 -n 5 -d 0 -T "${THREADS}" -M 0
cd-hit -i "${C99}"           -o "${C97}" -c 0.97 -aS 0.8 -g 1 -d 0 -T "${THREADS}" -M 0
cd-hit -i "${C97}"           -o "${C95}" -c 0.95 -aS 0.8 -g 1 -d 0 -T "${THREADS}" -M 0
cd-hit -i "${C95}"           -o "${C90}" -c 0.90 -aS 0.8 -g 1 -d 0 -T "${THREADS}" -M 0

# -------------------------
# 7. MAFFT alignment & Treemmer from 95% clustered set
# -------------------------
echo "[INFO] Running MAFFT on 95% clustered set..."
ALN95="${OUTDIR}/seq_c95.aln.fasta"
mafft --auto "${C95}" > "${ALN95}"

TREE95="${OUTDIR}/seq_c95.nwk"
echo "[INFO] Running FastTree (WAG+Gamma) on 95% alignment..."
FastTree -wag -gamma "${ALN95}" > "${TREE95}"

echo "[INFO] Running Treemmer to trim to N=${TARGET_N}..."
python "${TREEMMER_PY}" "${TREE95}" -X "${TARGET_N}" -np

# Treemmer v0.3 default naming
KEPT_LIST="${OUTDIR}/seq_c95.nwk_trimmed_list_X_${TARGET_N}"

if [[ ! -f "${KEPT_LIST}" ]]; then
    echo "[WARN] Treemmer kept-list not found at: ${KEPT_LIST}" >&2
    echo "      Please check Treemmer output naming; adjust script if needed." >&2
else
    echo "[INFO] Kept-list: ${KEPT_LIST}"
fi

REDUCED_ALN="${OUTDIR}/seq_${TARGET_N}_reduced.aln.fasta"
if [[ -f "${KEPT_LIST}" ]]; then
    seqkit grep -f "${KEPT_LIST}" "${ALN95}" > "${REDUCED_ALN}"
else
    echo "[WARN] Using full 95% alignment for final tree (no trimming)."
    cp "${ALN95}" "${REDUCED_ALN}"
fi

# -------------------------
# 8. Final FastTree
# -------------------------
FINAL_TREE="${OUTDIR}/seq_${TARGET_N}.tree"
echo "[INFO] Building final FastTree on reduced alignment..."
FastTree -wag -gamma "${REDUCED_ALN}" > "${FINAL_TREE}"

FINAL_IDS="${OUTDIR}/seq_${TARGET_N}_ids.txt"
grep '^>' "${REDUCED_ALN}" | sed 's/^>//' > "${FINAL_IDS}"

# -------------------------
# 9. Simple report
# -------------------------
REPORT="${OUTDIR}/report.txt"
{
    echo "=== Processing complete ==="
    echo "Input 1 fasta: ${INPUT1}"
    if [[ -n "${INPUT2}" ]]; then
        echo "Input 2 fasta: ${INPUT2}"
    fi
    echo "Min length (aa): ${MINLEN}"
    if [[ -n "${AA_FILTER}" ]]; then
        echo "Residue filter:  sequences containing at least one of [${AA_FILTER}]"
    else
        echo "Residue filter:  (none)"
    fi
    echo "After filters:    $(grep -c '^>' "${FILTERED_FASTA}") seqs"
    echo "After CD-HIT 99%: $(grep -c '^>' "${C99}") seqs"
    echo "After CD-HIT 97%: $(grep -c '^>' "${C97}") seqs"
    echo "After CD-HIT 95%: $(grep -c '^>' "${C95}") seqs"
    echo "After CD-HIT 90%: $(grep -c '^>' "${C90}") seqs"
    echo "Treemmer target N: ${TARGET_N}"
    echo "Final IDs:         ${FINAL_IDS}"
    echo "Final tree:        ${FINAL_TREE}"
} > "${REPORT}"

echo "[INFO] All done. See ${REPORT}"

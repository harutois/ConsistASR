#!/usr/bin/env bash
# ============================================================
# run_indel_aware_paml.sh  (RAxML + PAML pipeline)
#
# Indel-aware ASR pipeline for:
#   - RAxML tree used for PAML ASR (.nwk)
#   - PAML ASR output (rst)
#
# This script:
#   1) Converts an amino-acid MSA to a binary (0/1) gap matrix.
#   2) Evaluates the given topology with RAxML-NG under BIN+G
#      and writes a bestTree on the binary alignment.
#   3) Runs RAxML-HPC ancestral reconstruction on the binary
#      gap-state alignment.
#   4) Maps RAxML node IDs to PAML node IDs using the rst tree.
#   5) Writes raw PAML ASR FASTA from the rst file and merges
#      PAML rst AA states and indel states to produce gap-aware
#      ancestral FASTA.
#
# Example:
#   bash run_indel_aware_paml.sh \
#     --msa   HeR_SzR_228_FFT_NS_1.fasta \
#     --tree  HeR_SzR_228_FFT_NS_1.raxml.bestTree.tre \
#     --rst   HeR_SzR_228_FFT_NS_1.rst \
#     --prefix ASR_LGFG4_FFT_NS_1 \
#     --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
#
# Requirements:
#  - msa_to_binary.py, map_raxml_to_paml_nodes_from_rst.py, and
#    paml_state_and_indel_to_fasta.py in the same directory as this wrapper
#  - raxml-ng, raxmlHPC, and python in PATH
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MSA_FASTA=""
TREEFILE=""
RST_FILE=""
PREFIX=""
OUTGROUPS=""

usage() {
  cat <<EOF
Usage: $0 --msa MSA_FASTA --tree TREEFILE --rst RST_FILE --prefix OUT_PREFIX --outgroup "TAXA,..."

Required arguments:
  --msa       Amino-acid MSA in FASTA format (the exact alignment used for RAxML/PAML ASR)
  --tree      RAxML topology file (.nwk) used for the PAML ASR run
  --rst       PAML ASR output file (rst)
  --prefix    Output prefix (used for all intermediate and final files)
  --outgroup  Comma-separated list of outgroup taxa (passed to RAxML-NG --outgroup)

Example:
  $0 \\
    --msa   HeR_SzR_228_FFT_NS_1.fasta \\
    --tree  HeR_SzR_228_FFT_NS_1.raxml.bestTree.tre \\
    --rst   HeR_SzR_228_FFT_NS_1.rst \\
    --prefix ASR_LGFG4_FFT_NS_1 \\
    --outgroup "OG1,OG2,OG3"
EOF
}

if [ $# -eq 0 ]; then
   usage
   exit 1
fi

normalize_text_file() {
  local label="$1"
  local infile="$2"
  local outfile="$3"

  if [[ ! -f "$infile" ]]; then
    echo "[ERROR] $label file not found: $infile" >&2
    exit 1
  fi

  if grep -q $'\r' "$infile"; then
    echo "[WARN] CR/CRLF line endings detected in $label: $infile" >&2
    echo "[WARN] Creating LF-normalized working copy: $outfile" >&2
    tr -d '\r' < "$infile" > "$outfile"
  else
    cp "$infile" "$outfile"
  fi
}

# ---------------------------
# Parse command-line options
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --msa)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --msa requires an argument."
        usage
        exit 1
      fi
      MSA_FASTA="$2"
      shift 2
      ;;
    --tree)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --tree requires an argument."
        usage
        exit 1
      fi
      TREEFILE="$2"
      shift 2
      ;;
    --rst)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --rst requires an argument."
        usage
        exit 1
      fi
      RST_FILE="$2"
      shift 2
      ;;
    --prefix)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --prefix requires an argument."
        usage
        exit 1
      fi
      PREFIX="$2"
      shift 2
      ;;
    --outgroup)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --outgroup requires an argument."
        usage
        exit 1
      fi
      OUTGROUPS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# ---------------------------
# Basic argument checks
# ---------------------------
if [[ -z "$MSA_FASTA" || -z "$TREEFILE" || -z "$RST_FILE" || -z "$PREFIX" || -z "$OUTGROUPS" ]]; then
  echo "[ERROR] Missing required argument(s)."
  usage
  exit 1
fi

for f in "$MSA_FASTA" "$TREEFILE" "$RST_FILE"; do
  if [ ! -f "$f" ]; then
    echo "[ERROR] File not found: $f"
    exit 1
  fi
done

command -v raxml-ng   >/dev/null 2>&1 || { echo "[ERROR] raxml-ng not found in PATH";   exit 1; }
command -v raxmlHPC   >/dev/null 2>&1 || { echo "[ERROR] raxmlHPC not found in PATH";   exit 1; }
command -v python     >/dev/null 2>&1 || { echo "[ERROR] python not found in PATH";     exit 1; }

# These Python scripts are expected to be in the same directory as this wrapper:
#   - msa_to_binary.py
#   - map_raxml_to_paml_nodes_from_rst.py
#   - paml_state_and_indel_to_fasta.py

BIN_PHY="${PREFIX}_binary.phy"
EVAL_PREFIX="${PREFIX}_indel_eval"
RAxML_PREFIX="${PREFIX}_binary"
INDEL_PAML_NAMED="${PREFIX}_indel_ASR.paml_named.txt"
OUT_RAW_WITHGAP="${PREFIX}_raw_asr_withgap.fasta"
OUT_RAW_NOGAP="${PREFIX}_raw_asr_nogap.fasta"
OUT_NOGAP="${PREFIX}_indel_nogap.fasta"
OUT_WITHGAP="${PREFIX}_indel_withgap.fasta"
MAP_TABLE="${PREFIX}_node_map.tsv"
PAML_TREE="${PREFIX}_paml_nodes.tree"
WORK_SUBDIR="${PREFIX}_indel_work"

mkdir -p "$WORK_SUBDIR"

# ---------------------------
# Normalize line endings
# ---------------------------
# Original input files are not modified.
# LF-normalized working copies are used internally to avoid hidden
# Windows-style carriage returns in FASTA IDs, tree files, or rst files.

MSA_FASTA_ORIG="$MSA_FASTA"
TREEFILE_ORIG="$TREEFILE"
RST_FILE_ORIG="$RST_FILE"

MSA_FASTA_NORM="${WORK_SUBDIR}/input_msa.lf.fasta"
TREEFILE_NORM="${WORK_SUBDIR}/input_tree.lf.nwk"
RST_FILE_NORM="${WORK_SUBDIR}/input_paml_rst.lf"

normalize_text_file "MSA FASTA" "$MSA_FASTA_ORIG" "$MSA_FASTA_NORM"
normalize_text_file "RAxML tree" "$TREEFILE_ORIG" "$TREEFILE_NORM"
normalize_text_file "PAML rst" "$RST_FILE_ORIG" "$RST_FILE_NORM"

MSA_FASTA="$MSA_FASTA_NORM"
TREEFILE="$TREEFILE_NORM"
RST_FILE="$RST_FILE_NORM"

# Remove accidental CR/LF characters and spaces from comma-separated outgroup list.
OUTGROUPS="$(printf '%s' "$OUTGROUPS" | tr -d '\r\n' | sed 's/[[:space:]]//g')"

if [[ -z "$OUTGROUPS" ]]; then
  echo "[ERROR] --outgroup is empty after cleanup." >&2
  exit 1
fi

IFS=',' read -r -a OG_ARR <<< "$OUTGROUPS"

FASTA_IDS=$(grep '^>' "$MSA_FASTA" | sed 's/^>//' | awk '{print $1}')

for og in "${OG_ARR[@]}"; do
  if ! echo "$FASTA_IDS" | grep -qx "$og"; then
    echo "[ERROR] Outgroup taxon '$og' not found in MSA FASTA (>$og)" >&2
    echo "[HINT] Check FASTA sequence IDs and --outgroup list." >&2
    exit 1
  fi
done

echo "[INFO] =================================================="
echo "[INFO] 1) Converting MSA to binary (0/1): ${BIN_PHY}"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/msa_to_binary.py" "$MSA_FASTA" "$BIN_PHY"

echo "[INFO] =================================================="
echo "[INFO] 2) RAxML-NG --evaluate (BIN+G, outgroup=${OUTGROUPS})"
echo "[INFO] =================================================="
raxml-ng \
  --evaluate \
  --msa "$BIN_PHY" \
  --msa-format PHYLIP \
  --model BIN+G \
  --tree "$TREEFILE" \
  --prefix "$EVAL_PREFIX" \
  --outgroup "$OUTGROUPS" \
  --seed 12345

BEST_TREE="${EVAL_PREFIX}.raxml.bestTree"
if [ ! -f "$BEST_TREE" ]; then
  echo "[ERROR] Expected bestTree not found: $BEST_TREE"
  exit 1
fi

echo "[INFO] =================================================="
echo "[INFO] 3) RAxML-HPC -f A (BINGAMMA, binary gap-state ASR)"
echo "[INFO] =================================================="
raxmlHPC -f A \
  -m BINGAMMA \
  -t "$BEST_TREE" \
  -s "$BIN_PHY" \
  -n "$RAxML_PREFIX" \
  -F

RAxML_NODE_TREE="RAxML_nodeLabelledRootedTree.${RAxML_PREFIX}"
RAxML_ANCESTRAL="RAxML_marginalAncestralStates.${RAxML_PREFIX}"

if [ ! -f "$RAxML_NODE_TREE" ] || [ ! -f "$RAxML_ANCESTRAL" ]; then
  echo "[ERROR] Expected RAxML outputs not found:"
  echo "       $RAxML_NODE_TREE"
  echo "       $RAxML_ANCESTRAL"
  exit 1
fi

echo "[INFO] =================================================="
echo "[INFO] 4) Mapping node IDs (RAxML → PAML)"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/map_raxml_to_paml_nodes_from_rst.py" \
  --rst "$RST_FILE" \
  --raxml "$RAxML_NODE_TREE" \
  --ancestral "$RAxML_ANCESTRAL" \
  --out "$INDEL_PAML_NAMED" \
  --map_table "$MAP_TABLE" \
  --paml_tree_out "$PAML_TREE"

echo "[INFO] =================================================="
echo "[INFO] 5) Generating indel-aware ancestral FASTA"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/paml_state_and_indel_to_fasta.py" \
  --rst "$RST_FILE" \
  --indel "$INDEL_PAML_NAMED" \
  --out_raw_withgap "$OUT_RAW_WITHGAP" \
  --out_raw_nogap "$OUT_RAW_NOGAP" \
  --out_withgap "$OUT_WITHGAP" \
  --out_nogap "$OUT_NOGAP"

# --- 6) Move intermediate files to a subdirectory ---
echo "[INFO] =================================================="
echo "[INFO] Moving intermediate files to: ${WORK_SUBDIR}"
mkdir -p "$WORK_SUBDIR"

# Binary alignment and its reduced version
for f in "$BIN_PHY" "${BIN_PHY}.reduced"; do
  [[ -f "$f" ]] && mv -f "$f" "$WORK_SUBDIR"/
done

# raxml-ng intermediate files (evaluation on binary MSA)
mv -f ${EVAL_PREFIX}.raxml.* "$WORK_SUBDIR"/ 2>/dev/null || true

# RAxML-HPC outputs (info, tree, ancestral states, etc.)
mv -f RAxML_*."$RAxML_PREFIX" "$WORK_SUBDIR"/ 2>/dev/null || true

# Mapped binary indel states used for generating the final indel-aware FASTA
mv -f "$INDEL_PAML_NAMED" "$WORK_SUBDIR"/ 2>/dev/null || true

echo "[INFO] =================================================="
echo "[INFO] All done."
echo "[INFO]  - Indel-aware FASTA (with gaps)   : $OUT_WITHGAP"
echo "[INFO]  - Indel-aware FASTA (gap-stripped): $OUT_NOGAP"
echo "[INFO]  - raw FASTA (with gaps)           : $OUT_RAW_WITHGAP"
echo "[INFO]  - raw FASTA (gap-stripped)        : $OUT_RAW_NOGAP"
if [ -n "${MAP_TABLE:-}" ]; then
  echo "[INFO]  - Node mapping table              : $MAP_TABLE"
fi
echo "[INFO]  - PAML node-labelled tree         : $PAML_TREE"
echo "[INFO]  - Intermediate files moved to     : ${WORK_SUBDIR}"
echo "[INFO] =================================================="

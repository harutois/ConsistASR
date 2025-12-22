#!/usr/bin/env bash
# ============================================================
# run_indel_aware_paml.sh  (RAxML + PAML pipeline)
#
# Indel-aware ASR pipeline for:
#   - RAxML tree used for PAML ASR (.nwk)
#   - PAML ASR output (.rst)
#
# This script:
#   1) Converts an amino-acid MSA to a binary (0/1) gap matrix.
#   2) Evaluates the given topology with RAxML-NG under BIN+G
#      and writes a bestTree on the binary alignment.
#   3) Runs RAxML(-HPC) ancestral reconstruction on the binary
#      alignment (indel ASR).
#   4) Maps RAxML node IDs to PAML node IDs using the .rst tree.
#   5) Merges PAML .rst AA states and indel states to produce
#      gap-aware ancestral FASTA (with-gap and gap-stripped),
#      for all nodes.
#
# Example:
#   bash run_indel_aware_paml.sh \
#     --msa   HeR_SzR_228_FFT_NS_1.fasta \
#     --tree  HeR_SzR_228_FFT_NS_1.raxml.bestTree.tre \
#     --rst   HeR_SzR_228_FFT_NS_1.rst \
#     --prefix ASR_LGFG4_FFT_NS_1 \
#     --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
#
# Requirements (in PATH or same directory):
#   - msa_to_binary.py
#   - map_raxml_to_paml_nodes_from_rst.py
#   - paml_state_and_indel_to_fasta.py
#   - raxml-ng
#   - raxmlHPC
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
  --msa       Amino-acid MSA in FASTA format (same sequences as used for RAxML/PAML ASR)
  --tree      RAxML topology file (.nwk) used for the PAML ASR run
  --rst       PAML ASR output file (.rst)
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
  exit 1
}

if [ $# -eq 0 ]; then
  usage
fi

# ---------------------------
# Parse command-line options
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --msa)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --msa requires an argument."
        usage
      fi
      MSA_FASTA="$2"
      shift 2
      ;;
    --tree)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --tree requires an argument."
        usage
      fi
      TREEFILE="$2"
      shift 2
      ;;
    --rst)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --rst requires an argument."
        usage
      fi
      RST_FILE="$2"
      shift 2
      ;;
    --prefix)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --prefix requires an argument."
        usage
      fi
      PREFIX="$2"
      shift 2
      ;;
    --outgroup)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --outgroup requires an argument."
        usage
      fi
      OUTGROUPS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "[ERROR] Unknown option: $1"
      usage
      ;;
  esac
done

# ---------------------------
# Basic argument checks
# ---------------------------
if [[ -z "$MSA_FASTA" || -z "$TREEFILE" || -z "$RST_FILE" || -z "$PREFIX" || -z "$OUTGROUPS" ]]; then
  echo "[ERROR] Missing required argument(s)."
  usage
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

# These Python scripts are expected to be in the same directory or in PATH:
#   - msa_to_binary.py
#   - map_raxml_to_paml_nodes_from_rst.py
#   - paml_state_and_indel_to_fasta.py

BIN_PHY="${PREFIX}_binary.phy"
EVAL_PREFIX="${PREFIX}_indel_eval"
RAxML_PREFIX="${PREFIX}_binary"
INDEL_PAML_NAMED="${PREFIX}_indel_ASR.paml_named.txt"
OUT_NOGAP="${PREFIX}_indel_nogap.fasta"
OUT_WITHGAP="${PREFIX}_indel_withgap.fasta"
MAP_TABLE="${PREFIX}_node_map.tsv"
PAML_TREE="${PREFIX}_paml_nodes.tree"

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
echo "[INFO] 3) RAxML-HPC -f A (BINGAMMA, indel ASR)"
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
echo "[INFO] 5) Generating gap-aware ancestral FASTA (PAML ASR + indel)"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/paml_state_and_indel_to_fasta.py" \
  --rst "$RST_FILE" \
  --indel "$INDEL_PAML_NAMED" \
  --out_withgap "$OUT_WITHGAP" \
  --out_nogap "$OUT_NOGAP"

# --- 6) Move intermediate files to a subdirectory ---
WORK_SUBDIR="${PREFIX}_indel_work"
echo "[INFO] =================================================="
echo "[INFO] Moving intermediate files to: ${WORK_SUBDIR}"
mkdir -p "$WORK_SUBDIR"

# Binary alignment and its reduced version
mv -f "$BIN_PHY" "${BIN_PHY}.reduced" "$WORK_SUBDIR"/ 2>/dev/null || true

# raxml-ng intermediate files (evaluation on binary MSA)
mv -f ${EVAL_PREFIX}.raxml.* "$WORK_SUBDIR"/ 2>/dev/null || true

# RAxML-HPC outputs (info, tree, ancestral states, etc.)
mv -f RAxML_*."$RAxML_PREFIX" "$WORK_SUBDIR"/ 2>/dev/null || true

echo "[INFO] =================================================="
echo "[INFO] All done."
echo "[INFO]  - Indel-aware FASTA (with gaps):    $OUT_WITHGAP"
echo "[INFO]  - Indel-aware FASTA (gap-stripped): $OUT_NOGAP"
echo "[INFO]  - PAML node-labelled tree:          $PAML_TREE"
if [ -n "${MAP_TABLE:-}" ]; then
  echo "[INFO]  - Node mapping table:              $MAP_TABLE"
fi
echo "[INFO]  - Intermediate files moved to:     ${WORK_SUBDIR}"
echo "[INFO] =================================================="

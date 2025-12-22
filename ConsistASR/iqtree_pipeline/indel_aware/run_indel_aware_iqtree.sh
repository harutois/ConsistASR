#!/usr/bin/env bash
# ============================================================
# run_indel_aware_iqtree.sh
# Indel-aware ASR pipeline for IQ-TREE outputs (.state, .treefile).
#
# This script:
#   1) Converts an amino-acid MSA to a binary (0/1) gap matrix.
#   2) Evaluates the IQ-TREE topology with RAxML-NG under BIN+G
#      and writes a bestTree on the binary alignment.
#   3) Runs RAxML(-HPC) ancestral reconstruction on the binary
#      alignment (indel ASR).
#   4) Maps RAxML node IDs to IQ-TREE node IDs.
#   5) Merges IQ-TREE .state and indel states to produce
#      gap-aware ancestral FASTA (with-gap and gap-stripped).
#
# Example:
#   bash run_indel_aware_iqtree.sh \
#     --msa HeR_SzR_228_PSI_TM_Coffee.fasta \
#     --tree ASR_QpfamR7_PSI_TM_Coffee.treefile \
#     --state ASR_QpfamR7_PSI_TM_Coffee.state \
#     --prefix ASR_QpfamR7_PSI_TM_Coffee \
#     --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MSA_FASTA=""
TREEFILE=""
STATE_FILE=""
PREFIX=""
OUTGROUPS=""

usage() {
  cat <<EOF
Usage: $0 --msa MSA_FASTA --tree TREEFILE --state STATE_FILE --prefix OUT_PREFIX --outgroup "TAXA,..."

Required arguments:
  --msa       Amino-acid MSA in FASTA format (same sequences as used for IQ-TREE ASR)
  --tree      IQ-TREE topology file (.treefile) for the ASR run
  --state     IQ-TREE ancestral state file (.state)
  --prefix    Output prefix (used for all intermediate and final files)
  --outgroup  Comma-separated list of outgroup taxa (passed to RAxML-NG --outgroup)

Example:
  $0 \\
    --msa HeR_SzR_228_PSI_TM_Coffee.fasta \\
    --tree ASR_QpfamR7_PSI_TM_Coffee.treefile \\
    --state ASR_QpfamR7_PSI_TM_Coffee.state \\
    --prefix ASR_QpfamR7_PSI_TM_Coffee \\
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
    --state)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] --state requires an argument."
        usage
      fi
      STATE_FILE="$2"
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
if [[ -z "$MSA_FASTA" || -z "$TREEFILE" || -z "$STATE_FILE" || -z "$PREFIX" || -z "$OUTGROUPS" ]]; then
  echo "[ERROR] Missing required argument(s)."
  usage
fi

for f in "$MSA_FASTA" "$TREEFILE" "$STATE_FILE"; do
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
#   - map_raxml_to_iqtree_nodes.py
#   - state_and_indel_to_fasta.py

BIN_PHY="${PREFIX}_binary.phy"
EVAL_PREFIX="${PREFIX}_indel_eval"
RAxML_PREFIX="${PREFIX}_binary"
INDEL_MAP_OUT="${PREFIX}_indel_ASR.iqtree_named.txt"
OUT_NOGAP="${PREFIX}_indel_nogap.fasta"
OUT_WITHGAP="${PREFIX}_indel_withgap.fasta"
MAP_TABLE="${PREFIX}_node_map.tsv"

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
echo "[INFO] 4) Mapping node IDs (RAxML → IQ-TREE)"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/map_raxml_to_iqtree_nodes.py" \
  --iqtree "$TREEFILE" \
  --raxml "$RAxML_NODE_TREE" \
  --ancestral "$RAxML_ANCESTRAL" \
  --out "$INDEL_MAP_OUT" \
  --map_table "$MAP_TABLE"

echo "[INFO] =================================================="
echo "[INFO] 5) Generating gap-aware ancestral FASTA"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/state_and_indel_to_fasta.py" \
  --state "$STATE_FILE" \
  --indel "$INDEL_MAP_OUT" \
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
if [ -n "${MAP_TABLE:-}" ]; then
  echo "[INFO]  - Node mapping table:     　　　　$MAP_TABLE"
fi
echo "[INFO]  - Intermediate files moved to: 　　${WORK_SUBDIR}"
echo "[INFO] =================================================="

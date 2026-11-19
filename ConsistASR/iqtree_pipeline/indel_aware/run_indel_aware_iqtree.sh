#!/usr/bin/env bash
# ============================================================
# run_indel_aware_iqtree.sh
# Indel-aware ASR pipeline for IQ-TREE outputs (.state, .treefile).
#
# This script:
#   1) Converts an amino-acid MSA to a binary (0/1) gap matrix.
#   2) Evaluates the fixed IQ-TREE topology with RAxML-NG under BIN+G,
#      re-optimizing branch lengths and free model parameters for the
#      binary matrix.
#   3) Runs marginal binary ancestral-state reconstruction with
#      RAxML-NG --ancestral using the evaluated bestTree/bestModel.
#   4) Validates NodeXX/site correspondence against the IQ-TREE .state
#      file, derives binary calls from p_0/p_1, and writes a per-site
#      binary posterior-probability table.
#   5) Writes raw IQ-TREE ASR FASTA and merges IQ-TREE amino-acid states
#      with binary presence/absence states to produce indel-aware FASTA.
#
# Example:
#   bash run_indel_aware_iqtree.sh \
#     --msa HeR_SzR_228_PSI_TM_Coffee.fasta \
#     --tree ASR_QpfamR7_PSI_TM_Coffee.treefile \
#     --state ASR_QpfamR7_PSI_TM_Coffee.state \
#     --prefix ASR_QpfamR7_PSI_TM_Coffee \
#     --outgroup "OG_WP_285271495,OG_WP_136361479,OG_WP_010903286"
#
# Requirements:
#  - shared helpers msa_to_binary.py and raxmlng_binary_asr.py are either:
#      * in ConsistASR/common/, or
#      * copied alongside this wrapper for standalone use
#  - pipeline-specific helper state_and_indel_to_fasta.py is alongside this wrapper
#  - raxml-ng and python in PATH
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONSISTASR_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
COMMON_DIR="${CONSISTASR_ROOT}/common"

# Prefer colocated helpers when the workflow is copied as a standalone set.
# Otherwise, use the shared helpers in the repository-level common/ directory.
if [[ -f "${SCRIPT_DIR}/msa_to_binary.py" &&
      -f "${SCRIPT_DIR}/raxmlng_binary_asr.py" ]]; then
  SHARED_HELPER_DIR="${SCRIPT_DIR}"
  echo "[INFO] Using colocated shared helpers: ${SHARED_HELPER_DIR}"
elif [[ -f "${COMMON_DIR}/msa_to_binary.py" &&
        -f "${COMMON_DIR}/raxmlng_binary_asr.py" ]]; then
  SHARED_HELPER_DIR="${COMMON_DIR}"
  echo "[INFO] Using repository common helpers: ${SHARED_HELPER_DIR}"
else
  echo "[ERROR] Shared helper scripts were not found." >&2
  echo "[ERROR] Expected both msa_to_binary.py and raxmlng_binary_asr.py in either:" >&2
  echo "[ERROR]   ${SCRIPT_DIR}" >&2
  echo "[ERROR]   ${COMMON_DIR}" >&2
  exit 1
fi

MSA_FASTA=""
TREEFILE=""
STATE_FILE=""
PREFIX=""
OUTGROUPS=""

usage() {
  cat <<EOF
Usage: $0 --msa MSA_FASTA --tree TREEFILE --state STATE_FILE --prefix OUT_PREFIX --outgroup "TAXA,..."

Required arguments:
  --msa       Amino-acid MSA in FASTA format (the exact alignment used for IQ-TREE ASR)
  --tree      IQ-TREE topology file (.treefile) for the ASR run
  --state     IQ-TREE ancestral state file (.state)
  --prefix    Output prefix (used for all intermediate and final files)
  --outgroup  Comma-separated list of outgroup taxa (passed to RAxML-NG --outgroup)

Example:
  $0 \
    --msa HeR_SzR_228_PSI_TM_Coffee.fasta \
    --tree ASR_QpfamR7_PSI_TM_Coffee.treefile \
    --state ASR_QpfamR7_PSI_TM_Coffee.state \
    --prefix ASR_QpfamR7_PSI_TM_Coffee \
    --outgroup "OG1,OG2,OG3"
EOF
}

if [[ $# -eq 0 ]]; then
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
      [[ $# -ge 2 ]] || { echo "[ERROR] --msa requires an argument." >&2; usage; exit 1; }
      MSA_FASTA="$2"
      shift 2
      ;;
    --tree)
      [[ $# -ge 2 ]] || { echo "[ERROR] --tree requires an argument." >&2; usage; exit 1; }
      TREEFILE="$2"
      shift 2
      ;;
    --state)
      [[ $# -ge 2 ]] || { echo "[ERROR] --state requires an argument." >&2; usage; exit 1; }
      STATE_FILE="$2"
      shift 2
      ;;
    --prefix)
      [[ $# -ge 2 ]] || { echo "[ERROR] --prefix requires an argument." >&2; usage; exit 1; }
      PREFIX="$2"
      shift 2
      ;;
    --outgroup)
      [[ $# -ge 2 ]] || { echo "[ERROR] --outgroup requires an argument." >&2; usage; exit 1; }
      OUTGROUPS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# ---------------------------
# Basic argument checks
# ---------------------------
if [[ -z "$MSA_FASTA" || -z "$TREEFILE" || -z "$STATE_FILE" || -z "$PREFIX" || -z "$OUTGROUPS" ]]; then
  echo "[ERROR] Missing required argument(s)." >&2
  usage
  exit 1
fi

for f in "$MSA_FASTA" "$TREEFILE" "$STATE_FILE"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] File not found: $f" >&2
    exit 1
  fi
done

command -v raxml-ng >/dev/null 2>&1 || { echo "[ERROR] raxml-ng not found in PATH" >&2; exit 1; }
command -v python   >/dev/null 2>&1 || { echo "[ERROR] python not found in PATH" >&2; exit 1; }

if [[ ! -f "${SCRIPT_DIR}/state_and_indel_to_fasta.py" ]]; then
  echo "[ERROR] Required pipeline-specific script not found: ${SCRIPT_DIR}/state_and_indel_to_fasta.py" >&2
  exit 1
fi

BIN_PHY="${PREFIX}_binary.phy"
EVAL_PREFIX="${PREFIX}_indel_eval"
ASR_PREFIX="${PREFIX}_indel_asr"

BINARY_PROFILE="${PREFIX}_binary_states.txt"
BINARY_PP_TSV="${PREFIX}_binary_posterior.tsv"

OUT_RAW_WITHGAP="${PREFIX}_raw_asr_withgap.fasta"
OUT_RAW_NOGAP="${PREFIX}_raw_asr_nogap.fasta"
OUT_NOGAP="${PREFIX}_indel_nogap.fasta"
OUT_WITHGAP="${PREFIX}_indel_withgap.fasta"
OUT_SITE_POSTERIOR="${PREFIX}_site_posterior.tsv"
OUT_NOGAP_POSTERIOR="${PREFIX}_indel_nogap_posterior.tsv"

WORK_SUBDIR="${PREFIX}_indel_work"
mkdir -p "$WORK_SUBDIR"

# ---------------------------
# Normalize line endings
# ---------------------------
# Original input files are not modified.
MSA_FASTA_ORIG="$MSA_FASTA"
TREEFILE_ORIG="$TREEFILE"
STATE_FILE_ORIG="$STATE_FILE"

MSA_FASTA_NORM="${WORK_SUBDIR}/input_msa.lf.fasta"
TREEFILE_NORM="${WORK_SUBDIR}/input_tree.lf.nwk"
STATE_FILE_NORM="${WORK_SUBDIR}/input_iqtree_state.lf.state"

normalize_text_file "MSA FASTA" "$MSA_FASTA_ORIG" "$MSA_FASTA_NORM"
normalize_text_file "IQ-TREE tree" "$TREEFILE_ORIG" "$TREEFILE_NORM"
normalize_text_file "IQ-TREE .state" "$STATE_FILE_ORIG" "$STATE_FILE_NORM"

MSA_FASTA="$MSA_FASTA_NORM"
TREEFILE="$TREEFILE_NORM"
STATE_FILE="$STATE_FILE_NORM"

# Clean comma-separated outgroup list.
OUTGROUPS="$(printf '%s' "$OUTGROUPS" | tr -d '\r\n' | sed 's/[[:space:]]//g')"

if [[ -z "$OUTGROUPS" ]]; then
  echo "[ERROR] --outgroup is empty after cleanup." >&2
  exit 1
fi

IFS=',' read -r -a OG_ARR <<< "$OUTGROUPS"
FASTA_IDS="$(grep '^>' "$MSA_FASTA" | sed 's/^>//' | awk '{print $1}')"

for og in "${OG_ARR[@]}"; do
  if ! printf '%s\n' "$FASTA_IDS" | grep -qx "$og"; then
    echo "[ERROR] Outgroup taxon '$og' not found in MSA FASTA (>$og)" >&2
    echo "[HINT] Check FASTA sequence IDs and --outgroup list." >&2
    exit 1
  fi
done

echo "[INFO] =================================================="
echo "[INFO] 1) Converting MSA to binary (0/1): ${BIN_PHY}"
echo "[INFO] =================================================="
python "${SHARED_HELPER_DIR}/msa_to_binary.py" "$MSA_FASTA" "$BIN_PHY"

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
BEST_MODEL="${EVAL_PREFIX}.raxml.bestModel"

for f in "$BEST_TREE" "$BEST_MODEL"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Expected RAxML-NG evaluation output not found: $f" >&2
    exit 1
  fi
done

echo "[INFO] =================================================="
echo "[INFO] 3) RAxML-NG --ancestral (evaluated BIN+G model)"
echo "[INFO] =================================================="
raxml-ng \
  --ancestral \
  --msa "$BIN_PHY" \
  --msa-format PHYLIP \
  --model "$BEST_MODEL" \
  --tree "$BEST_TREE" \
  --prefix "$ASR_PREFIX" \
  --outgroup "$OUTGROUPS" \
  --opt-model off \
  --opt-branches off \
  --seed 12345

BINARY_PROBS="${ASR_PREFIX}.raxml.ancestralProbs"
BINARY_STATES="${ASR_PREFIX}.raxml.ancestralStates"

for f in "$BINARY_PROBS" "$BINARY_STATES"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Expected RAxML-NG ancestral output not found: $f" >&2
    exit 1
  fi
done

echo "[INFO] =================================================="
echo "[INFO] 4) Validating binary ASR and writing posterior table"
echo "[INFO] =================================================="
python "${SHARED_HELPER_DIR}/raxmlng_binary_asr.py" \
  --probs "$BINARY_PROBS" \
  --states "$BINARY_STATES" \
  --iqtree-state "$STATE_FILE" \
  --out-binary "$BINARY_PROFILE" \
  --out-tsv "$BINARY_PP_TSV"

echo "[INFO] =================================================="
echo "[INFO] 5) Generating indel-aware ancestral FASTA"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/state_and_indel_to_fasta.py" \
  --state "$STATE_FILE" \
  --indel "$BINARY_PROFILE" \
  --binary_posterior "$BINARY_PP_TSV" \
  --out_raw_withgap "$OUT_RAW_WITHGAP" \
  --out_raw_nogap "$OUT_RAW_NOGAP" \
  --out_withgap "$OUT_WITHGAP" \
  --out_nogap "$OUT_NOGAP" \
  --out_site_posterior "$OUT_SITE_POSTERIOR" \
  --out_nogap_posterior "$OUT_NOGAP_POSTERIOR"

echo "[INFO] =================================================="
echo "[INFO] 6) Moving intermediate files to: ${WORK_SUBDIR}"
echo "[INFO] =================================================="

# Binary alignment and any RAxML-NG reduced copy.
for f in "$BIN_PHY" "${BIN_PHY}.reduced"; do
  [[ -f "$f" ]] && mv -f "$f" "$WORK_SUBDIR"/
done

# RAxML-NG evaluation outputs.
for f in ${EVAL_PREFIX}.raxml.*; do
  [[ -e "$f" ]] && mv -f "$f" "$WORK_SUBDIR"/
done

# RAxML-NG ancestral outputs.
for f in ${ASR_PREFIX}.raxml.*; do
  [[ -e "$f" ]] && mv -f "$f" "$WORK_SUBDIR"/
done

# Internal two-column binary profiles used for FASTA merging.
[[ -f "$BINARY_PROFILE" ]] && mv -f "$BINARY_PROFILE" "$WORK_SUBDIR"/

echo "[INFO] =================================================="
echo "[INFO] All done."
echo "[INFO]  - Indel-aware FASTA (with gaps)   : $OUT_WITHGAP"
echo "[INFO]  - Indel-aware FASTA (gap-stripped): $OUT_NOGAP"
echo "[INFO]  - raw FASTA (with gaps)           : $OUT_RAW_WITHGAP"
echo "[INFO]  - raw FASTA (gap-stripped)        : $OUT_RAW_NOGAP"
echo "[INFO]  - Binary posterior table          : $BINARY_PP_TSV"
echo "[INFO]  - Site posterior table            : $OUT_SITE_POSTERIOR"
echo "[INFO]  - Final nogap posterior table     : $OUT_NOGAP_POSTERIOR"
echo "[INFO]  - Intermediate files moved to     : ${WORK_SUBDIR}"
echo "[INFO] =================================================="

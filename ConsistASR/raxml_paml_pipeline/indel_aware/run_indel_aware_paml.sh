#!/usr/bin/env bash
# ============================================================
# run_indel_aware_paml.sh  (RAxML-NG + PAML pipeline)
#
# Indel-aware ASR pipeline for:
#   - the amino-acid tree used for PAML ASR
#   - PAML ASR output (rst)
#
# The PAML degree-2 root has no one-to-one RAxML-NG binary-ASR
# counterpart and is therefore excluded from indel-aware outputs.
# It remains available in raw PAML FASTA and the complete site-level
# posterior table.
#
# Steps:
#   1) Convert the amino-acid MSA to a binary 0/1 gap matrix.
#   2) RAxML-NG --evaluate under BIN+G (topology retained; branch
#      lengths/model parameters optimized on binary data).
#   3) RAxML-NG --ancestral using the evaluated tree/model with
#      model/branch optimization disabled.
#   4) Validate binary ancestral probabilities and derive 0/1 calls.
#   5) Map all n-2 RAxML-NG internal vertices one-to-one to PAML
#      non-root ancestral node IDs; explicitly exclude the PAML root.
#   6) Merge PAML amino-acid states with mapped binary states and write
#      FASTA plus posterior-probability TSVs.
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
RST_FILE=""
PREFIX=""
OUTGROUPS=""

usage() {
  cat <<EOF
Usage: $0 --msa MSA_FASTA --tree TREEFILE --rst RST_FILE --prefix OUT_PREFIX --outgroup "TAXA,..."

Required arguments:
  --msa       Amino-acid MSA in FASTA format (exact alignment used for PAML ASR)
  --tree      Amino-acid tree/topology file used for the PAML ASR run
  --rst       PAML ASR output file (rst)
  --prefix    Output prefix
  --outgroup  Comma-separated outgroup taxa passed to RAxML-NG

Example:
  $0 \\
    --msa toy_7tm_rhodopsin.fasta \\
    --tree ASR_toy_7tm_rhodopsin_paml_input.nwk \\
    --rst rst \\
    --prefix ASR_raxmlng_paml_toy_7tm_rhodopsin \\
    --outgroup "OG_WP_285271495,OG_WP_136361479"
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

while [[ $# -gt 0 ]]; do
  case "$1" in
    --msa)
      [[ $# -ge 2 ]] || { echo "[ERROR] --msa requires an argument."; exit 1; }
      MSA_FASTA="$2"; shift 2 ;;
    --tree)
      [[ $# -ge 2 ]] || { echo "[ERROR] --tree requires an argument."; exit 1; }
      TREEFILE="$2"; shift 2 ;;
    --rst)
      [[ $# -ge 2 ]] || { echo "[ERROR] --rst requires an argument."; exit 1; }
      RST_FILE="$2"; shift 2 ;;
    --prefix)
      [[ $# -ge 2 ]] || { echo "[ERROR] --prefix requires an argument."; exit 1; }
      PREFIX="$2"; shift 2 ;;
    --outgroup)
      [[ $# -ge 2 ]] || { echo "[ERROR] --outgroup requires an argument."; exit 1; }
      OUTGROUPS="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      exit 1 ;;
  esac
done

if [[ -z "$MSA_FASTA" || -z "$TREEFILE" || -z "$RST_FILE" || -z "$PREFIX" || -z "$OUTGROUPS" ]]; then
  echo "[ERROR] Missing required argument(s)." >&2
  usage
  exit 1
fi

for f in "$MSA_FASTA" "$TREEFILE" "$RST_FILE"; do
  [[ -f "$f" ]] || { echo "[ERROR] File not found: $f" >&2; exit 1; }
done

command -v raxml-ng >/dev/null 2>&1 || { echo "[ERROR] raxml-ng not found in PATH" >&2; exit 1; }
command -v python   >/dev/null 2>&1 || { echo "[ERROR] python not found in PATH" >&2; exit 1; }

for s in map_raxml_to_paml_nodes_from_rst.py paml_state_and_indel_to_fasta.py; do
  [[ -f "${SCRIPT_DIR}/${s}" ]] || {
    echo "[ERROR] Required pipeline-specific script not found: ${SCRIPT_DIR}/${s}" >&2
    exit 1
  }
done

BIN_PHY="${PREFIX}_binary.phy"
EVAL_PREFIX="${PREFIX}_indel_eval"
ASR_PREFIX="${PREFIX}_indel_asr"

RAXML_BINARY="${PREFIX}_binary_states.raxmlng.txt"
RAXML_BINARY_PP="${PREFIX}_binary_posterior.raxmlng.tsv"

INDEL_PAML_NAMED="${PREFIX}_indel_ASR.paml_named.txt"
BINARY_PP_PAML="${PREFIX}_binary_posterior.tsv"

OUT_RAW_WITHGAP="${PREFIX}_raw_asr_withgap.fasta"
OUT_RAW_NOGAP="${PREFIX}_raw_asr_nogap.fasta"
OUT_WITHGAP="${PREFIX}_indel_withgap.fasta"
OUT_NOGAP="${PREFIX}_indel_nogap.fasta"

OUT_SITE_POSTERIOR="${PREFIX}_site_posterior.tsv"
OUT_NOGAP_POSTERIOR="${PREFIX}_indel_nogap_posterior.tsv"

MAP_TABLE="${PREFIX}_node_map.tsv"
PAML_TREE="${PREFIX}_paml_nodes.tree"
WORK_SUBDIR="${PREFIX}_indel_work"

mkdir -p "$WORK_SUBDIR"

MSA_FASTA_NORM="${WORK_SUBDIR}/input_msa.lf.fasta"
TREEFILE_NORM="${WORK_SUBDIR}/input_tree.lf.nwk"
RST_FILE_NORM="${WORK_SUBDIR}/input_paml_rst.lf"

normalize_text_file "MSA FASTA" "$MSA_FASTA" "$MSA_FASTA_NORM"
normalize_text_file "tree" "$TREEFILE" "$TREEFILE_NORM"
normalize_text_file "PAML rst" "$RST_FILE" "$RST_FILE_NORM"

MSA_FASTA="$MSA_FASTA_NORM"
TREEFILE="$TREEFILE_NORM"
RST_FILE="$RST_FILE_NORM"

OUTGROUPS="$(printf '%s' "$OUTGROUPS" | tr -d '\r\n' | sed 's/[[:space:]]//g')"
[[ -n "$OUTGROUPS" ]] || { echo "[ERROR] --outgroup is empty after cleanup." >&2; exit 1; }

IFS=',' read -r -a OG_ARR <<< "$OUTGROUPS"
FASTA_IDS=$(grep '^>' "$MSA_FASTA" | sed 's/^>//' | awk '{print $1}')

for og in "${OG_ARR[@]}"; do
  if ! echo "$FASTA_IDS" | grep -qx "$og"; then
    echo "[ERROR] Outgroup taxon '$og' not found in MSA FASTA (>$og)" >&2
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
  [[ -f "$f" ]] || { echo "[ERROR] Expected RAxML-NG output not found: $f" >&2; exit 1; }
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

ASR_PROBS="${ASR_PREFIX}.raxml.ancestralProbs"
ASR_STATES="${ASR_PREFIX}.raxml.ancestralStates"
ASR_TREE="${ASR_PREFIX}.raxml.ancestralTree"

for f in "$ASR_PROBS" "$ASR_STATES" "$ASR_TREE"; do
  [[ -f "$f" ]] || { echo "[ERROR] Expected RAxML-NG ancestral output not found: $f" >&2; exit 1; }
done

echo "[INFO] =================================================="
echo "[INFO] 4) Validating RAxML-NG binary ASR"
echo "[INFO] =================================================="
python "${SHARED_HELPER_DIR}/raxmlng_binary_asr.py" \
  --probs "$ASR_PROBS" \
  --states "$ASR_STATES" \
  --out-binary "$RAXML_BINARY" \
  --out-tsv "$RAXML_BINARY_PP"

echo "[INFO] =================================================="
echo "[INFO] 5) Mapping RAxML-NG nodes to PAML non-root nodes"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/map_raxml_to_paml_nodes_from_rst.py" \
  --rst "$RST_FILE" \
  --raxml-tree "$ASR_TREE" \
  --binary "$RAXML_BINARY" \
  --posterior "$RAXML_BINARY_PP" \
  --out-binary "$INDEL_PAML_NAMED" \
  --out-posterior "$BINARY_PP_PAML" \
  --map-table "$MAP_TABLE" \
  --paml-tree-out "$PAML_TREE"

PAML_ROOT="$(
  awk -F '\t' 'NR>1 && $3=="excluded_paml_root" {print $1}' "$MAP_TABLE"
)"

if [[ -z "$PAML_ROOT" ]]; then
  echo "[ERROR] PAML root node was not recorded in mapping table." >&2
  exit 1
fi

if [[ "$(printf '%s\n' "$PAML_ROOT" | wc -l | tr -d ' ')" -ne 1 ]]; then
  echo "[ERROR] Expected exactly one excluded PAML root node." >&2
  exit 1
fi

echo "[INFO] PAML root excluded from indel-aware integration: ${PAML_ROOT}"

echo "[INFO] =================================================="
echo "[INFO] 6) Generating PAML indel-aware FASTA and posterior tables"
echo "[INFO] =================================================="
python "${SCRIPT_DIR}/paml_state_and_indel_to_fasta.py" \
  --rst "$RST_FILE" \
  --indel "$INDEL_PAML_NAMED" \
  --binary-posterior "$BINARY_PP_PAML" \
  --root-node "$PAML_ROOT" \
  --out_raw_withgap "$OUT_RAW_WITHGAP" \
  --out_raw_nogap "$OUT_RAW_NOGAP" \
  --out_withgap "$OUT_WITHGAP" \
  --out_nogap "$OUT_NOGAP" \
  --out-site-posterior "$OUT_SITE_POSTERIOR" \
  --out-nogap-posterior "$OUT_NOGAP_POSTERIOR"

echo "[INFO] =================================================="
echo "[INFO] 7) Moving intermediate files to: ${WORK_SUBDIR}"
echo "[INFO] =================================================="

for f in "$BIN_PHY" "${BIN_PHY}.reduced" "$RAXML_BINARY" "$RAXML_BINARY_PP" "$INDEL_PAML_NAMED"; do
  [[ -f "$f" ]] && mv -f "$f" "$WORK_SUBDIR"/
done

mv -f ${EVAL_PREFIX}.raxml.* "$WORK_SUBDIR"/ 2>/dev/null || true
mv -f ${ASR_PREFIX}.raxml.* "$WORK_SUBDIR"/ 2>/dev/null || true

echo "[INFO] =================================================="
echo "[INFO] All done."
echo "[INFO]  - Indel-aware FASTA (with gaps)   : $OUT_WITHGAP"
echo "[INFO]  - Indel-aware FASTA (gap-stripped): $OUT_NOGAP"
echo "[INFO]  - raw PAML FASTA (with gaps)       : $OUT_RAW_WITHGAP"
echo "[INFO]  - raw PAML FASTA (gap-stripped)    : $OUT_RAW_NOGAP"
echo "[INFO]  - Binary posterior table (PAML IDs): $BINARY_PP_PAML"
echo "[INFO]  - Site posterior table             : $OUT_SITE_POSTERIOR"
echo "[INFO]  - Final nogap posterior table      : $OUT_NOGAP_POSTERIOR"
echo "[INFO]  - Node mapping table               : $MAP_TABLE"
echo "[INFO]  - PAML node-labelled tree          : $PAML_TREE"
echo "[INFO]  - Excluded PAML root node          : $PAML_ROOT"
echo "[INFO]  - Intermediate files moved to      : ${WORK_SUBDIR}"
echo "[INFO] =================================================="

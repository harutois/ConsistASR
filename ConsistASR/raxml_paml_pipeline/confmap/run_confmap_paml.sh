#!/usr/bin/env bash
# ============================================================
# run_confmap_paml.sh
#
# Wrapper script to combine site-wise posterior probabilities
# from a PAML .rst file with an indel-aware ancestral alignment
# (withgap.fasta) and AlphaFold models, and embed these values
# into the B-factors of a CIF/PDB structure.
#
# The script:
#   1. Automatically detects AlphaFold JSON/CIF files for the target node
#   2. Converts CIF → PDB (assuming pLDDT is stored in the B-factor)
#   3. Extracts per-alignment-column PP values for the specified node
#      from the PAML .rst file
#   4. Extracts the corresponding with-gap alignment sequence for that node
#      from the withgap FASTA
#      (no indel masking or expansion is performed here; handling of gap
#       columns is delegated to map_confidence_to_bfactor.py)
#   5. Uses map_confidence_to_bfactor.py to embed PP into B-factors
#   6. Generates additional PDB files where the B-factor stores:
#        - scaled PP − pLDDT
#        - PP × pLDDT
#   7. Computes summary statistics (mean, min, max, median) for PP, pLDDT,
#      PP − pLDDT, and PP × pLDDT (CA-only) and writes them to a log file
#
# Requirements:
#   - extract_pp_from_paml_rst.py
#   - map_confidence_to_bfactor.py
#   must be available in $PATH or in the same directory.
#
# Usage:
#   bash run_confmap_paml.sh \
#       --rst     <PAML_RST> \
#       --node    <PAML_NODE_INT> \
#       --withgap <INDEL_WITHGAP_FASTA> \
#       [--outdir confmap]
#
# Example:
#   bash run_confmap_paml.sh \
#       --rst     HeR_SzR_228_FFT_NS_1.rst \
#       --node    233 \
#       --withgap ASR_LGFG4_FFT_NS_1_indel_withgap.fasta
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# -------------------------
# Argument parsing
# -------------------------
RST_FILE=""
PAML_NODE=""
WITHGAP_FASTA=""
OUTDIR="confmap"

usage() {
    echo "Usage: $0 --rst <PAML_RST> --node <PAML_NODE_INT> --withgap <INDEL_WITHGAP_FASTA> [--outdir confmap]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --rst)
            RST_FILE="$2"; shift 2;;
        --node)
            PAML_NODE="$2"; shift 2;;
        --withgap)
            WITHGAP_FASTA="$2"; shift 2;;
        --outdir)
            OUTDIR="$2"; shift 2;;
        -h|--help)
            usage; exit 0;;
        *)
            echo "[ERROR] Unknown option: $1"
            usage
            exit 1;;
    esac
done

if [[ -z "${RST_FILE}" || -z "${PAML_NODE}" || -z "${WITHGAP_FASTA}" ]]; then
    echo "[ERROR] --rst, --node, --withgap は必須です。"
    usage
    exit 1
fi

if [[ ! -f "${RST_FILE}" ]]; then
    echo "[ERROR] RST_FILE not found: ${RST_FILE}"
    exit 1
fi

if [[ ! -f "${WITHGAP_FASTA}" ]]; then
    echo "[ERROR] WITHGAP_FASTA not found: ${WITHGAP_FASTA}"
    exit 1
fi

NODE_NUM="${PAML_NODE}"
NODE_LABEL="Node${NODE_NUM}"

export NODE_NAME="${NODE_LABEL}"

# ------------------------------------------------------------
# Create output directory (e.g., confmap/)
# (Assumes this script is run inside the AlphaFold NodeXXX dir)
# ------------------------------------------------------------
mkdir -p "${OUTDIR}"
OUTDIR_ABS="$(cd "${OUTDIR}" && pwd)"

echo "[INFO] Output directory: ${OUTDIR_ABS}"
echo "[INFO] PAML node: ${PAML_NODE} (label=${NODE_LABEL})"

# Define output file paths (all under outputdirectory)
PLDDT_PDB="${OUTDIR}/${NODE_LABEL}_plddt_bfactor.pdb"
TMP_PP_FILE="${OUTDIR}/${NODE_LABEL}_pp_values.txt"
NODE_ALN_FASTA="${OUTDIR}/${NODE_LABEL}_aligned_indel.fasta"
PP_PDB="${OUTDIR}/${NODE_LABEL}_pp_bfactor.pdb"
DIFF_PDB="${OUTDIR}/${NODE_LABEL}_ppminusplddt_bfactor.pdb"
PROD_PDB="${OUTDIR}/${NODE_LABEL}_ppxplddt_bfactor.pdb"
STATS_LOG="${OUTDIR}/${NODE_LABEL}_conf_CA_stats.log"

# ------------------------------------------------------------
# Step 1) Detect AlphaFold JSON/CIF for this node
# ------------------------------------------------------------
echo "[INFO] --- Step 1: Searching AlphaFold model outputs ---"

BEST_INFO=$(python - <<'PY'
import glob, json, re, os, sys

node = os.environ.get("NODE_NAME", "")
node_lower = node.lower()

candidates = sorted(glob.glob("*_summary_confidences_*.json"))
if not candidates:
    print("[ERROR] No *_summary_confidences_*.json found in current directory.", file=sys.stderr)
    sys.exit(1)

# Prefer files whose name contains the node id (if any)
node_files = [f for f in candidates if node_lower and node_lower in f.lower()]
files = node_files or candidates

def score_file(path):
    try:
        with open(path) as f:
            d = json.load(f)
        val = d.get("ranking_score", d.get("plddt", 0.0))
        return float(val)
    except Exception:
        return 0.0

best = max(files, key=score_file)
best_score = score_file(best)

m = re.search(r"_summary_confidences_(\d+)\.json$", best)
if m:
    idx = int(m.group(1))
else:
    idx = 0

print(best)
print(idx)
print(best_score)
PY
)

if [ -z "${BEST_INFO}" ]; then
    echo "[ERROR] Failed to select best summary_confidences JSON."
    exit 1
fi

BEST_JSON=$(echo "${BEST_INFO}" | sed -n '1p')
BEST_INDEX=$(echo "${BEST_INFO}" | sed -n '2p')
BEST_SCORE=$(echo "${BEST_INFO}" | sed -n '3p')

echo "[INFO] Detected JSON: ${BEST_JSON}"
echo "[INFO] Best model index (from filename): ${BEST_INDEX}"
echo "[INFO] Score value (ranking_score or pLDDT): ${BEST_SCORE}"

PREFIX="${BEST_JSON%_summary_confidences_*}"
BEST_CIF="${PREFIX}_model_${BEST_INDEX}.cif"

if [ ! -f "${BEST_CIF}" ]; then
    echo "[ERROR] CIF not found: ${BEST_CIF}"
    exit 1
fi
echo "[INFO] Detected CIF: ${BEST_CIF}"

# ------------------------------------------------------------
# Step 2) CIF → PDB conversion (keep pLDDT in B-factor)
# ------------------------------------------------------------
echo "[INFO] --- Step 2: Converting CIF to PDB with pLDDT ---"

if command -v phenix.cif_as_pdb &> /dev/null; then
    echo "[INFO] Using phenix.cif_as_pdb"
    phenix.cif_as_pdb "${BEST_CIF}" > "${PLDDT_PDB}"
elif command -v gemmi &> /dev/null; then
    echo "[INFO] Using gemmi"
    gemmi convert "${BEST_CIF}" "${PLDDT_PDB}"
elif command -v pdb_fromcif &> /dev/null; then
    echo "[INFO] Using pdb_fromcif"
    pdb_fromcif "${BEST_CIF}" > "${PLDDT_PDB}"
else
    echo "[ERROR] No CIF→PDB converter found (phenix / gemmi / pdb_fromcif)"
    exit 1
fi
echo "[INFO] pLDDT PDB: ${PLDDT_PDB}"

# ------------------------------------------------------------
# Step 3) Extract PP values for this node from rst
# ------------------------------------------------------------
echo "[INFO] --- Step 3: Extracting PP from PAML RST for node ${PAML_NODE} ---"

python "${SCRIPT_DIR}/extract_pp_from_paml_rst.py" \
  --rst "${RST_FILE}" \
  --paml_node "${PAML_NODE}" \
  --withgap_fasta "${WITHGAP_FASTA}" \
  --out "${TMP_PP_FILE}"

# ------------------------------------------------------------
# Step 4) Extract aligned sequence for this node (indel-aware)
# ------------------------------------------------------------
echo "[INFO] --- Step 4: Extracting aligned sequence for node ${PAML_NODE} from withgap FASTA ---"

awk -v node=">${PAML_NODE}" '
  $0==node { print $0; getline; print; exit }
' "${WITHGAP_FASTA}" > "${NODE_ALN_FASTA}"

if [ ! -s "${NODE_ALN_FASTA}" ]; then
    echo "[ERROR] >${PAML_NODE} not found in ${WITHGAP_FASTA} (exact header match)."
    exit 1
fi
echo "[INFO] Node alignment fasta: ${NODE_ALN_FASTA}"

# ------------------------------------------------------------
# Step 5) Embed PP into B-factors (PP in 0–100 scale)
# ------------------------------------------------------------
echo "[INFO] --- Step 5: Embedding PP into B-factors via map_confidence_to_bfactor.py ---"

python "${SCRIPT_DIR}/map_confidence_to_bfactor.py" \
  --cif "${BEST_CIF}" \
  --asr_fasta "${NODE_ALN_FASTA}" \
  --pp_file "${TMP_PP_FILE}" \
  --mode pp \
  --out "${PP_PDB}"

echo "[INFO] PP-mapped PDB: ${PP_PDB}"

# ------------------------------------------------------------
# Step 6) Generate PDBs for PP−pLDDT (scaled) and PP×pLDDT
# ------------------------------------------------------------
echo "[INFO] --- Step 6: Generating PP - pLDDT and PP x pLDDT PDBs ---"

python - <<PY
from Bio.PDB import PDBParser, PDBIO
from copy import deepcopy

plddt_pdb = "${PLDDT_PDB}"
pp_pdb    = "${PP_PDB}"
diff_pdb  = "${DIFF_PDB}"
prod_pdb  = "${PROD_PDB}"

parser = PDBParser(QUIET=True)
pl_structure = parser.get_structure("PL", plddt_pdb)
pp_structure = parser.get_structure("PP", pp_pdb)

pl_atoms = list(pl_structure.get_atoms())
pp_atoms = list(pp_structure.get_atoms())

n = min(len(pl_atoms), len(pp_atoms))
pl_atoms = pl_atoms[:n]
pp_atoms = pp_atoms[:n]

diff_structure = deepcopy(pp_structure)
diff_atoms = list(diff_structure.get_atoms())

prod_structure = deepcopy(pp_structure)
prod_atoms = list(prod_structure.get_atoms())

for i, (pa, pl) in enumerate(zip(pp_atoms, pl_atoms)):
    pp = pa.get_bfactor()
    plv = pl.get_bfactor()

    # diff: (PP - pLDDT) in [-100, 100] → scale to [0, 100]
    diff = pp - plv
    scaled_diff = (diff + 100.0) / 200.0 * 100.0
    diff_atoms[i].set_bfactor(scaled_diff)

    # prod: (PP * pLDDT / 100) in [0, 100]
    prod = (pp * plv) / 100.0
    prod_atoms[i].set_bfactor(prod)

io = PDBIO()
io.set_structure(diff_structure)
io.save(diff_pdb)

io.set_structure(prod_structure)
io.save(prod_pdb)

print(f"[INFO] Wrote diff PDB (PP - pLDDT scaled): {diff_pdb}")
print(f"[INFO] Wrote prod PDB (PP x pLDDT):       {prod_pdb}")
PY

# ------------------------------------------------------------
# Step 7) CA-only statistics for PP / pLDDT / diff / prod
# ------------------------------------------------------------
echo "[INFO] --- Step 7: Computing CA-only statistics for PP / pLDDT / diff / prod ---"

python - <<PY
import statistics
from Bio.PDB import PDBParser

node_label = "${NODE_LABEL}"
pp_file = "${TMP_PP_FILE}"
plddt_pdb = "${PLDDT_PDB}"
pp_pdb    = "${PP_PDB}"
diff_pdb  = "${DIFF_PDB}"
prod_pdb  = "${PROD_PDB}"
log_file  = "${STATS_LOG}"

def stats(arr):
    if not arr:
        return (0.0, 0.0, 0.0, 0.0)
    return (
        statistics.mean(arr),
        min(arr),
        max(arr),
        statistics.median(arr),
    )

# PP_sitewise from .state (residue-level, one per alignment column)
pp_vals = []
with open(pp_file) as f:
    for line in f:
        s = line.strip()
        if not s:
            continue
        try:
            pp_vals.append(float(s))
        except ValueError:
            pass
pp_mean, pp_min, pp_max, pp_med = stats(pp_vals)

parser = PDBParser(QUIET=True)
def bfacts_ca(path):
    vals = []
    structure = parser.get_structure("X", path)
    for atom in structure.get_atoms():
        if atom.get_name() == "CA":
            vals.append(atom.get_bfactor())
    return vals

pl_vals   = bfacts_ca(plddt_pdb)
pp_bvals  = bfacts_ca(pp_pdb)
diff_vals = bfacts_ca(diff_pdb)
prod_vals = bfacts_ca(prod_pdb)

pl_mean, pl_min, pl_max, pl_med     = stats(pl_vals)
ppb_mean, ppb_min, ppb_max, ppb_med = stats(pp_bvals)
df_mean, df_min, df_max, df_med     = stats(diff_vals)
pr_mean, pr_min, pr_max, pr_med     = stats(prod_vals)

with open(log_file, "w") as lf:
    lf.write(f"[INFO] Node {node_label} per-residue (CA-only) statistics\\n")
    lf.write(f"PP_sitewise: mean={pp_mean:.3f}, min={pp_min:.3f}, max={pp_max:.3f}, median={pp_med:.3f}, n={len(pp_vals)}\\n")
    lf.write(f"pLDDT_CA:    mean={pl_mean:.3f}, min={pl_min:.3f}, max={pl_max:.3f}, median={pl_med:.3f}, n={len(pl_vals)}\\n")
    lf.write(f"PP_CA:       mean={ppb_mean:.3f}, min={ppb_min:.3f}, max={ppb_max:.3f}, median={ppb_med:.3f}, n={len(pp_bvals)}\\n")
    lf.write(f"PP-PLDDT_CA: mean={df_mean:.3f}, min={df_min:.3f}, max={df_max:.3f}, median={df_med:.3f}, n={len(diff_vals)}\\n")
    lf.write(f"PPxPLDDT_CA: mean={pr_mean:.3f}, min={pr_min:.3f}, max={pr_max:.3f}, median={pr_med:.3f}, n={len(prod_vals)}\\n")

print(f"[INFO] Log saved: {log_file}")
PY

echo "[INFO] --- Done ---"
echo "[INFO] PDB with pLDDT:       ${PLDDT_PDB}"
echo "[INFO] PDB with PP:          ${PP_PDB}"
echo "[INFO] PDB with PP-pLDDT:    ${DIFF_PDB}"
echo "[INFO] PDB with PPxPLDDT:    ${PROD_PDB}"
echo "[INFO] Stats log:            ${STATS_LOG}"

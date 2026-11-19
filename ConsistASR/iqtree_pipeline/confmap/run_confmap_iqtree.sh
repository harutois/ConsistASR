#!/usr/bin/env bash
# ============================================================
# run_confmap_iqtree.sh
#
# Map IQ-TREE ASR posterior probabilities (PP) onto AlphaFold
# models as B-factors, and generate additional PDB files where
# B-factors represent PP×pLDDT and scaled (PP−pLDDT).
#
# Workflow:
#   1. Automatically detect the AlphaFold JSON/CIF pair for a
#      given internal node from *_summary_confidences_*.json.
#   2. Convert the selected CIF to PDB, keeping pLDDT in the
#      B-factor field.
#   3. Extract site-wise PP values for the specified node from
#      the IQ-TREE .state file (maximum posterior per site).
#      Site numbering is explicitly checked.
#   4. Extract the gap-containing aligned ancestral sequence
#      for this node from the supplied with-gap FASTA.
#      Multi-line FASTA entries are supported.
#   5. Validate alignment-column count, structure length, and
#      structure/ancestral-sequence identity.
#   6. Use map_confidence_to_bfactor.py to embed PP as
#      per-residue B-factors in a new PDB.
#   7. Compute PP−pLDDT (scaled to 0–100) and PP×pLDDT
#      (0–100), using residue identifiers rather than atom order.
#   8. Compute CA-only statistics and write a log file.
#
# Usage:
#   bash run_confmap_iqtree.sh \
#       --state     <STATE_FILE> \
#       --node      <NODE> \
#       --withgap   <WITHGAP_FASTA> \
#       [--outdir confmap]
#
# Example:
#   bash run_confmap_iqtree.sh \
#       --state     ASR_QpfamR7_PSI_TM_Coffee.state \
#       --node      Node112 \
#       --withgap   ASR_QpfamR7_PSI_TM_Coffee_indel_withgap.fasta
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONSISTASR_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
COMMON_DIR="${CONSISTASR_ROOT}/common"

# Prefer a colocated copy for standalone use; otherwise use common/.
if [[ -f "${SCRIPT_DIR}/map_confidence_to_bfactor.py" ]]; then
  MAP_SCRIPT="${SCRIPT_DIR}/map_confidence_to_bfactor.py"
elif [[ -f "${COMMON_DIR}/map_confidence_to_bfactor.py" ]]; then
  MAP_SCRIPT="${COMMON_DIR}/map_confidence_to_bfactor.py"
else
  echo "[ERROR] map_confidence_to_bfactor.py not found." >&2
  echo "[ERROR] Expected it in either:" >&2
  echo "[ERROR]   ${SCRIPT_DIR}" >&2
  echo "[ERROR]   ${COMMON_DIR}" >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"

# -------------------------
# Argument parsing
# -------------------------
STATE_FILE=""
NODE=""
WITHGAP_FASTA=""
OUTDIR="confmap"

usage() {
    cat <<EOF
Usage:
  $0 --state <STATE_FILE> --node <NODE> --withgap <WITHGAP_FASTA> [--outdir confmap]

Required:
  --state PATH      IQ-TREE .state file
  --node STR        Internal-node identifier (e.g. Node112)
  --withgap PATH    Gap-containing ancestral FASTA preserving the
                    alignment coordinates of the IQ-TREE .state file

Optional:
  --outdir DIR      Output directory (default: confmap)
  -h, --help        Show this help and exit
EOF
}

require_value() {
    local opt="$1"
    local val="${2-}"

    if [[ -z "${val}" || "${val}" == -* ]]; then
        echo "[ERROR] ${opt} requires an argument." >&2
        usage >&2
        exit 1
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --state)
            require_value "$1" "${2-}"
            STATE_FILE="$2"
            shift 2
            ;;
        --node)
            require_value "$1" "${2-}"
            NODE="$2"
            shift 2
            ;;
        --withgap)
            require_value "$1" "${2-}"
            WITHGAP_FASTA="$2"
            shift 2
            ;;
        --outdir)
            require_value "$1" "${2-}"
            OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${STATE_FILE}" || -z "${NODE}" || -z "${WITHGAP_FASTA}" ]]; then
    echo "[ERROR] --state, --node, and --withgap are required." >&2
    usage >&2
    exit 1
fi

# -------------------------
# Basic checks
# -------------------------
if [[ ! -f "${STATE_FILE}" ]]; then
    echo "[ERROR] STATE_FILE not found: ${STATE_FILE}" >&2
    exit 1
fi

if [[ ! -f "${WITHGAP_FASTA}" ]]; then
    echo "[ERROR] WITHGAP_FASTA not found: ${WITHGAP_FASTA}" >&2
    exit 1
fi

if [[ ! -f "${MAP_SCRIPT}" ]]; then
    echo "[ERROR] Required helper script not found: ${MAP_SCRIPT}" >&2
    exit 1
fi

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    echo "[ERROR] Python executable not found: ${PYTHON_BIN}" >&2
    exit 1
fi

export NODE_NAME="${NODE}"

# ------------------------------------------------------------
# Create output directory
# (Assumes this script is run inside the AlphaFold node directory)
# ------------------------------------------------------------
mkdir -p "${OUTDIR}"
OUTDIR_ABS="$(cd "${OUTDIR}" && pwd)"

echo "[INFO] Output directory: ${OUTDIR_ABS}"

# Define output file paths
PLDDT_PDB="${OUTDIR}/${NODE}_plddt_bfactor.pdb"
TMP_PP_FILE="${OUTDIR}/${NODE}_pp_values.txt"
NODE_ALN_FASTA="${OUTDIR}/${NODE}_aligned_withgap.fasta"
PP_PDB="${OUTDIR}/${NODE}_pp_bfactor.pdb"
DIFF_PDB="${OUTDIR}/${NODE}_ppminusplddt_bfactor.pdb"
PROD_PDB="${OUTDIR}/${NODE}_ppxplddt_bfactor.pdb"
STATS_LOG="${OUTDIR}/${NODE}_conf_CA_stats.log"

# ------------------------------------------------------------
# Step 1) Detect AlphaFold JSON/CIF for this node
# ------------------------------------------------------------
echo "[INFO] --- Step 1: Searching AlphaFold model outputs ---"

BEST_INFO=$("${PYTHON_BIN}" - <<'PY'
import glob
import json
import math
import os
import re
import sys

# ------------------------------------------------------------
# Collect AlphaFold Server summary-confidence JSON files.
#
# The user-defined AlphaFold job name is intentionally NOT used
# to identify the ancestral node. The wrapper assumes that the
# current directory contains the output of one AlphaFold run
# corresponding to the sequence specified by --node/--withgap.
# ------------------------------------------------------------

candidates = sorted(glob.glob("*_summary_confidences_*.json"))

if not candidates:
    sys.exit(
        "[ERROR] No *_summary_confidences_*.json found "
        "in the current directory."
    )

# ------------------------------------------------------------
# Parse each filename into:
#   prefix + model index
#
# Example:
# fold_..._my_job_summary_confidences_0.json
# ->
# prefix = fold_..._my_job
# index  = 0
# ------------------------------------------------------------

parsed = []

for path in candidates:
    m = re.match(
        r"^(.*)_summary_confidences_(\d+)\.json$",
        os.path.basename(path),
    )

    if not m:
        print(
            f"[WARN] Skipping unrecognized summary filename: {path}",
            file=sys.stderr,
        )
        continue

    prefix = m.group(1)
    index = int(m.group(2))

    parsed.append((path, prefix, index))

if not parsed:
    sys.exit(
        "[ERROR] No summary-confidence filename could be parsed."
    )

# ------------------------------------------------------------
# Require exactly one AlphaFold job prefix.
#
# This avoids silently mixing outputs from different AlphaFold
# runs placed in the same directory.
# ------------------------------------------------------------

prefixes = sorted({prefix for _, prefix, _ in parsed})

if len(prefixes) != 1:
    msg = [
        "[ERROR] Multiple AlphaFold output prefixes were found "
        "in the current directory.",
        "The confidence-mapping wrapper expects outputs from "
        "one AlphaFold job per working directory.",
        "Detected prefixes:",
    ]

    msg.extend(f"  - {prefix}" for prefix in prefixes)

    sys.exit("\n".join(msg))

run_prefix = prefixes[0]

# ------------------------------------------------------------
# Validate JSON/CIF pairs and extract ranking scores.
# ------------------------------------------------------------

valid = []

for json_path, prefix, index in parsed:
    cif_path = f"{prefix}_model_{index}.cif"

    if not os.path.isfile(cif_path):
        print(
            f"[WARN] Skipping JSON without corresponding CIF: "
            f"{json_path} -> {cif_path}",
            file=sys.stderr,
        )
        continue

    try:
        with open(json_path) as fh:
            data = json.load(fh)
    except Exception as exc:
        print(
            f"[WARN] Skipping unreadable JSON {json_path}: {exc}",
            file=sys.stderr,
        )
        continue

    score_value = data.get("ranking_score")

    if score_value is None:
        score_value = data.get("plddt")

    if score_value is None:
        print(
            f"[WARN] Skipping JSON without ranking_score or plddt: "
            f"{json_path}",
            file=sys.stderr,
        )
        continue

    try:
        score = float(score_value)
    except (TypeError, ValueError):
        print(
            f"[WARN] Skipping JSON with non-numeric score: "
            f"{json_path}",
            file=sys.stderr,
        )
        continue

    if not math.isfinite(score):
        print(
            f"[WARN] Skipping JSON with non-finite score: "
            f"{json_path}",
            file=sys.stderr,
        )
        continue

    valid.append(
        (score, json_path, cif_path, index)
    )

if not valid:
    sys.exit(
        "[ERROR] No valid AlphaFold summary-confidence JSON/CIF "
        "pair was found."
    )

# Highest ranking score wins.
best_score, best_json, best_cif, best_index = max(
    valid,
    key=lambda x: x[0],
)

print(best_json)
print(best_cif)
print(best_index)
print(best_score)
PY
)

if [[ -z "${BEST_INFO}" ]]; then
    echo "[ERROR] Failed to select best summary-confidence JSON." >&2
    exit 1
fi

BEST_JSON=$(printf '%s\n' "${BEST_INFO}" | sed -n '1p')
BEST_CIF=$(printf '%s\n' "${BEST_INFO}" | sed -n '2p')
BEST_INDEX=$(printf '%s\n' "${BEST_INFO}" | sed -n '3p')
BEST_SCORE=$(printf '%s\n' "${BEST_INFO}" | sed -n '4p')

echo "[INFO] Detected JSON: ${BEST_JSON}"
echo "[INFO] Detected CIF:  ${BEST_CIF}"
echo "[INFO] Best model index: ${BEST_INDEX}"
echo "[INFO] Score value (ranking_score or pLDDT): ${BEST_SCORE}"

echo "[INFO] Detected CIF: ${BEST_CIF}"

# ------------------------------------------------------------
# Step 2) CIF -> PDB conversion (keep pLDDT in B-factor)
# ------------------------------------------------------------
echo "[INFO] --- Step 2: Converting CIF to PDB with pLDDT ---"

if command -v phenix.cif_as_pdb >/dev/null 2>&1; then
    echo "[INFO] Using phenix.cif_as_pdb"
    phenix.cif_as_pdb "${BEST_CIF}" > "${PLDDT_PDB}"

elif command -v gemmi >/dev/null 2>&1; then
    echo "[INFO] Using gemmi"
    gemmi convert "${BEST_CIF}" "${PLDDT_PDB}"

elif command -v pdb_fromcif >/dev/null 2>&1; then
    echo "[INFO] Using pdb_fromcif"
    pdb_fromcif "${BEST_CIF}" > "${PLDDT_PDB}"

else
    echo "[ERROR] No CIF->PDB converter found (phenix / gemmi / pdb_fromcif)." >&2
    exit 1
fi

if [[ ! -s "${PLDDT_PDB}" ]]; then
    echo "[ERROR] CIF-to-PDB conversion produced an empty file: ${PLDDT_PDB}" >&2
    exit 1
fi

echo "[INFO] pLDDT PDB: ${PLDDT_PDB}"

# ------------------------------------------------------------
# Step 3) Extract PP values for this node from .state
# ------------------------------------------------------------
echo "[INFO] --- Step 3: Extracting PP values from .state for ${NODE} ---"

"${PYTHON_BIN}" - \
    "${STATE_FILE}" \
    "${NODE}" \
    "${TMP_PP_FILE}" <<'PY'
import sys

state_file, node, out_file = sys.argv[1:]

pp_by_site = {}

with open(state_file) as fh:
    for lineno, raw in enumerate(fh, start=1):
        line = raw.strip()

        if not line or line.startswith("#"):
            continue

        parts = line.split()

        if not parts or parts[0] != node:
            continue

        if len(parts) < 4:
            sys.exit(
                f"[ERROR] Malformed .state row for {node} "
                f"at line {lineno}: expected at least 4 columns."
            )

        try:
            site = int(parts[1])
        except ValueError:
            sys.exit(
                f"[ERROR] Invalid site number for {node} "
                f"at line {lineno}: {parts[1]!r}"
            )

        if site <= 0:
            sys.exit(
                f"[ERROR] Invalid non-positive site number for "
                f"{node}: {site}"
            )

        if site in pp_by_site:
            sys.exit(
                f"[ERROR] Duplicate .state site {site} found "
                f"for node {node}."
            )

        try:
            vals = [float(x) for x in parts[3:]]
        except ValueError:
            sys.exit(
                f"[ERROR] Non-numeric posterior probability "
                f"for {node} at site {site}."
            )

        if not vals:
            sys.exit(
                f"[ERROR] No posterior probabilities found "
                f"for {node} at site {site}."
            )

        pp_by_site[site] = max(vals)

if not pp_by_site:
    sys.exit(
        f"[ERROR] No PP rows found for node {node} in {state_file}"
    )

sites = sorted(pp_by_site)

expected = list(range(1, sites[-1] + 1))

if sites != expected:
    missing = sorted(set(expected) - set(sites))
    preview = ", ".join(map(str, missing[:10]))

    if len(missing) > 10:
        preview += ", ..."

    sys.exit(
        f"[ERROR] Non-contiguous .state site numbering for {node}. "
        f"Missing site(s): {preview}"
    )

with open(out_file, "w") as out:
    for site in sites:
        out.write(f"{pp_by_site[site]:.6f}\n")

print(
    f"[INFO] Extracted {len(sites)} PP values for {node} "
    f"(sites 1-{sites[-1]}) -> {out_file}"
)
PY

# ------------------------------------------------------------
# Step 4) Extract aligned sequence for this node
# ------------------------------------------------------------
echo "[INFO] --- Step 4: Extracting aligned sequence for ${NODE} ---"

"${PYTHON_BIN}" - \
    "${WITHGAP_FASTA}" \
    "${NODE}" \
    "${NODE_ALN_FASTA}" <<'PY'
import sys

fasta, node, out_file = sys.argv[1:]

matches = []
header = None
seq_parts = []


def store_current():
    if header == node:
        matches.append("".join(seq_parts))


with open(fasta) as fh:
    for lineno, raw in enumerate(fh, start=1):
        line = raw.rstrip("\r\n")

        if not line.strip():
            continue

        if line.startswith(">"):
            store_current()
            header = line[1:]
            seq_parts = []
        else:
            if header is None:
                sys.exit(
                    f"[ERROR] Sequence data found before the first "
                    f"FASTA header at line {lineno} in {fasta}."
                )

            seq_parts.append("".join(line.split()))

store_current()

if len(matches) == 0:
    sys.exit(
        f"[ERROR] Exact FASTA header >{node} not found in {fasta}."
    )

if len(matches) > 1:
    sys.exit(
        f"[ERROR] Duplicate exact FASTA header >{node} found in {fasta}."
    )

sequence = matches[0]

if not sequence:
    sys.exit(
        f"[ERROR] Empty sequence found for >{node} in {fasta}."
    )

with open(out_file, "w") as out:
    out.write(f">{node}\n")
    out.write(sequence + "\n")

print(
    f"[INFO] Extracted aligned sequence for {node}: "
    f"{len(sequence)} alignment columns -> {out_file}"
)
PY

echo "[INFO] Node alignment FASTA: ${NODE_ALN_FASTA}"

# ------------------------------------------------------------
# Step 5) Validate PP/alignment/structure correspondence
# ------------------------------------------------------------
echo "[INFO] --- Step 5: Validating sequence and coordinate correspondence ---"

"${PYTHON_BIN}" - \
    "${NODE_ALN_FASTA}" \
    "${TMP_PP_FILE}" \
    "${PLDDT_PDB}" \
    "${NODE}" <<'PY'
import sys

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

fasta_file, pp_file, pdb_file, node = sys.argv[1:]

# -------------------------
# Read aligned ancestral sequence
# -------------------------
seq_parts = []

with open(fasta_file) as fh:
    for raw in fh:
        line = raw.strip()

        if not line or line.startswith(">"):
            continue

        seq_parts.append(line)

aligned = "".join(seq_parts).upper()

if not aligned:
    sys.exit("[ERROR] Extracted ancestral sequence is empty.")

ungapped = aligned.replace("-", "")

# -------------------------
# Read PP values
# -------------------------
pp_vals = []

with open(pp_file) as fh:
    for lineno, raw in enumerate(fh, start=1):
        s = raw.strip()

        if not s:
            continue

        try:
            pp_vals.append(float(s))
        except ValueError:
            sys.exit(
                f"[ERROR] Invalid PP value at line {lineno} "
                f"in {pp_file}: {s!r}"
            )

if len(pp_vals) != len(aligned):
    sys.exit(
        f"[ERROR] Alignment-column mismatch for {node}: "
        f"{len(pp_vals)} PP sites in .state versus "
        f"{len(aligned)} columns in the with-gap ancestral sequence."
    )

# -------------------------
# Extract protein sequence from AlphaFold PDB
# -------------------------
parser = PDBParser(QUIET=True)
structure = parser.get_structure("AF", pdb_file)

models = list(structure)

if len(models) != 1:
    sys.exit(
        f"[ERROR] Expected a single AlphaFold model in {pdb_file}, "
        f"but found {len(models)} models."
    )

protein_chains = []

for chain in models[0]:
    residues = []

    for residue in chain:
        if not is_aa(residue, standard=False):
            continue

        if "CA" not in residue:
            sys.exit(
                f"[ERROR] Protein residue {residue.get_resname()} "
                f"{residue.id[1]} in chain {chain.id!r} "
                f"has no CA atom."
            )

        aa = seq1(
            residue.get_resname(),
            undef_code="X",
        )

        residues.append(aa)

    if residues:
        protein_chains.append(
            (chain.id, "".join(residues).upper())
        )

if not protein_chains:
    sys.exit(
        f"[ERROR] No protein chain containing CA atoms "
        f"was found in {pdb_file}."
    )

if len(protein_chains) != 1:
    desc = ", ".join(
        f"{chain_id!r}:{len(seq)} aa"
        for chain_id, seq in protein_chains
    )

    sys.exit(
        f"[ERROR] Expected a monomeric AlphaFold structure, "
        f"but found {len(protein_chains)} protein chains ({desc})."
    )

chain_id, structure_seq = protein_chains[0]

if len(ungapped) != len(structure_seq):
    sys.exit(
        f"[ERROR] Ungapped ancestral length != structure sequence length "
        f"for {node}: ancestral={len(ungapped)}, "
        f"structure={len(structure_seq)}."
    )

if ungapped != structure_seq:
    mismatch = None

    for i, (a, b) in enumerate(
        zip(ungapped, structure_seq),
        start=1,
    ):
        if a != b:
            mismatch = (i, a, b)
            break

    if mismatch is None:
        detail = "sequence mismatch detected"
    else:
        pos, anc_aa, struct_aa = mismatch
        detail = (
            f"first mismatch at residue {pos}: "
            f"ancestral={anc_aa}, structure={struct_aa}"
        )

    sys.exit(
        f"[ERROR] AlphaFold structure sequence does not exactly match "
        f"the ungapped ancestral sequence for {node}; {detail}."
    )

print(
    f"[INFO] PP/alignment columns: {len(pp_vals)} / {len(aligned)}"
)
print(
    f"[INFO] Ungapped ancestral sequence: {len(ungapped)} residues"
)
print(
    f"[INFO] AlphaFold structure chain {chain_id!r}: "
    f"{len(structure_seq)} residues"
)
print(
    "[INFO] Structure sequence exactly matches the "
    "ungapped ancestral sequence."
)
PY

# ------------------------------------------------------------
# Step 6) Embed PP into B-factors (PP in 0-100 scale)
# ------------------------------------------------------------
echo "[INFO] --- Step 6: Embedding PP into B-factors ---"

"${PYTHON_BIN}" "${MAP_SCRIPT}" \
    --cif "${BEST_CIF}" \
    --asr_fasta "${NODE_ALN_FASTA}" \
    --pp_file "${TMP_PP_FILE}" \
    --mode pp \
    --out "${PP_PDB}"

if [[ ! -s "${PP_PDB}" ]]; then
    echo "[ERROR] PP-mapped PDB was not generated: ${PP_PDB}" >&2
    exit 1
fi

echo "[INFO] PP-mapped PDB: ${PP_PDB}"

# ------------------------------------------------------------
# Step 7) Generate PDBs for PP-pLDDT and PPxpLDDT
#
# IMPORTANT:
# Residues are matched by model / chain / residue identifier /
# residue name. Atom-list order is not used for correspondence.
# ------------------------------------------------------------
echo "[INFO] --- Step 7: Generating PP-pLDDT and PPxpLDDT PDBs ---"

"${PYTHON_BIN}" - \
    "${PLDDT_PDB}" \
    "${PP_PDB}" \
    "${DIFF_PDB}" \
    "${PROD_PDB}" <<'PY'
import sys
from copy import deepcopy

from Bio.PDB import PDBIO, PDBParser

plddt_pdb, pp_pdb, diff_pdb, prod_pdb = sys.argv[1:]

parser = PDBParser(QUIET=True)

pl_structure = parser.get_structure("PL", plddt_pdb)
pp_structure = parser.get_structure("PP", pp_pdb)


def residue_key(model, chain, residue):
    hetflag, resseq, icode = residue.id

    return (
        model.id,
        chain.id,
        hetflag,
        resseq,
        icode,
        residue.get_resname(),
    )


def ca_bfactor_map(structure, label):
    values = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" not in residue:
                    continue

                key = residue_key(model, chain, residue)

                if key in values:
                    raise RuntimeError(
                        f"Duplicate residue key in {label}: {key}"
                    )

                values[key] = residue["CA"].get_bfactor()

    if not values:
        raise RuntimeError(
            f"No CA-containing residues found in {label}."
        )

    return values


try:
    pl_map = ca_bfactor_map(pl_structure, "pLDDT PDB")
    pp_map = ca_bfactor_map(pp_structure, "PP PDB")
except RuntimeError as exc:
    sys.exit(f"[ERROR] {exc}")

pl_keys = set(pl_map)
pp_keys = set(pp_map)

if pl_keys != pp_keys:
    only_pl = sorted(pl_keys - pp_keys, key=str)
    only_pp = sorted(pp_keys - pl_keys, key=str)

    msg = [
        "[ERROR] Residue identifiers differ between "
        "the pLDDT and PP PDB files."
    ]

    if only_pl:
        msg.append(
            "Present only in pLDDT PDB: "
            + ", ".join(map(str, only_pl[:5]))
        )

    if only_pp:
        msg.append(
            "Present only in PP PDB: "
            + ", ".join(map(str, only_pp[:5]))
        )

    sys.exit("\n".join(msg))

# Validate confidence ranges before computing derived metrics.
for key in sorted(pp_keys, key=str):
    pp = pp_map[key]
    pl = pl_map[key]

    if not (0.0 <= pp <= 100.0):
        sys.exit(
            f"[ERROR] PP B-factor outside 0-100 range "
            f"at residue {key}: {pp}"
        )

    if not (0.0 <= pl <= 100.0):
        sys.exit(
            f"[ERROR] pLDDT B-factor outside 0-100 range "
            f"at residue {key}: {pl}"
        )

diff_values = {}
prod_values = {}

for key in pp_keys:
    pp = pp_map[key]
    pl = pl_map[key]

    # Raw difference: [-100, +100]
    diff = pp - pl

    # Scaled difference: [0, 100]
    scaled_diff = (diff + 100.0) / 2.0

    # Product metric: [0, 100]
    prod = (pp * pl) / 100.0

    diff_values[key] = scaled_diff
    prod_values[key] = prod


def apply_values(structure, values):
    seen = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" not in residue:
                    continue

                key = residue_key(model, chain, residue)

                if key not in values:
                    raise RuntimeError(
                        f"Residue missing from metric map: {key}"
                    )

                value = values[key]

                for atom in residue:
                    atom.set_bfactor(value)

                seen.add(key)

    if seen != set(values):
        missing = set(values) - seen

        raise RuntimeError(
            f"{len(missing)} residue(s) were not written "
            "to the output structure."
        )


diff_structure = deepcopy(pp_structure)
prod_structure = deepcopy(pp_structure)

try:
    apply_values(diff_structure, diff_values)
    apply_values(prod_structure, prod_values)
except RuntimeError as exc:
    sys.exit(f"[ERROR] {exc}")

io = PDBIO()

io.set_structure(diff_structure)
io.save(diff_pdb)

io.set_structure(prod_structure)
io.save(prod_pdb)

print(
    f"[INFO] Matched {len(pp_keys)} residues "
    "between PP and pLDDT structures."
)
print(
    f"[INFO] Wrote diff PDB (PP - pLDDT scaled): {diff_pdb}"
)
print(
    f"[INFO] Wrote prod PDB (PP x pLDDT):       {prod_pdb}"
)
PY

# ------------------------------------------------------------
# Step 8) CA-only statistics
# ------------------------------------------------------------
echo "[INFO] --- Step 8: Computing CA-only statistics ---"

"${PYTHON_BIN}" - \
    "${TMP_PP_FILE}" \
    "${PLDDT_PDB}" \
    "${PP_PDB}" \
    "${DIFF_PDB}" \
    "${PROD_PDB}" \
    "${STATS_LOG}" \
    "${NODE}" <<'PY'
import statistics
import sys

from Bio.PDB import PDBParser

(
    pp_file,
    plddt_pdb,
    pp_pdb,
    diff_pdb,
    prod_pdb,
    log_file,
    node,
) = sys.argv[1:]


def stats(arr, label):
    if not arr:
        sys.exit(
            f"[ERROR] No values available for statistics: {label}"
        )

    return (
        statistics.mean(arr),
        min(arr),
        max(arr),
        statistics.median(arr),
    )


# -------------------------
# PP_sitewise from .state
# -------------------------
pp_vals_raw = []

with open(pp_file) as fh:
    for lineno, raw in enumerate(fh, start=1):
        s = raw.strip()

        if not s:
            continue

        try:
            pp_vals_raw.append(float(s))
        except ValueError:
            sys.exit(
                f"[ERROR] Invalid PP value at line {lineno} "
                f"in {pp_file}: {s!r}"
            )

if not pp_vals_raw:
    sys.exit("[ERROR] No site-wise PP values were loaded.")

if min(pp_vals_raw) < 0:
    sys.exit("[ERROR] Negative PP value detected.")

# IQ-TREE normally reports probabilities on the 0-1 scale.
# Retain compatibility with inputs already expressed as 0-100.
if max(pp_vals_raw) <= 1.0:
    pp_vals = [v * 100.0 for v in pp_vals_raw]
elif max(pp_vals_raw) <= 100.0:
    pp_vals = pp_vals_raw
else:
    sys.exit(
        "[ERROR] Site-wise PP values exceed the expected 0-100 range."
    )

pp_mean, pp_min, pp_max, pp_med = stats(
    pp_vals,
    "PP_sitewise",
)

# -------------------------
# CA B-factors
# -------------------------
parser = PDBParser(QUIET=True)


def bfacts_ca(path):
    structure = parser.get_structure(path, path)

    vals = [
        atom.get_bfactor()
        for atom in structure.get_atoms()
        if atom.get_name() == "CA"
    ]

    return vals


pl_vals = bfacts_ca(plddt_pdb)
pp_bvals = bfacts_ca(pp_pdb)
diff_vals = bfacts_ca(diff_pdb)
prod_vals = bfacts_ca(prod_pdb)

n_ca = {
    "pLDDT": len(pl_vals),
    "PP": len(pp_bvals),
    "PP-pLDDT": len(diff_vals),
    "PPxpLDDT": len(prod_vals),
}

if len(set(n_ca.values())) != 1:
    sys.exit(
        "[ERROR] CA counts differ among confidence-mapped PDB files: "
        + ", ".join(
            f"{name}={count}"
            for name, count in n_ca.items()
        )
    )

pl_mean, pl_min, pl_max, pl_med = stats(
    pl_vals,
    "pLDDT_CA",
)

ppb_mean, ppb_min, ppb_max, ppb_med = stats(
    pp_bvals,
    "PP_CA",
)

df_mean, df_min, df_max, df_med = stats(
    diff_vals,
    "PP-pLDDT_CA",
)

pr_mean, pr_min, pr_max, pr_med = stats(
    prod_vals,
    "PPxpLDDT_CA",
)

with open(log_file, "w") as lf:
    lf.write(
        f"[INFO] Node {node} confidence statistics\n"
    )

    lf.write(
        "PP_sitewise_0to100: "
        f"mean={pp_mean:.3f}, "
        f"min={pp_min:.3f}, "
        f"max={pp_max:.3f}, "
        f"median={pp_med:.3f}, "
        f"n={len(pp_vals)}\n"
    )

    lf.write(
        "pLDDT_CA:            "
        f"mean={pl_mean:.3f}, "
        f"min={pl_min:.3f}, "
        f"max={pl_max:.3f}, "
        f"median={pl_med:.3f}, "
        f"n={len(pl_vals)}\n"
    )

    lf.write(
        "PP_CA:               "
        f"mean={ppb_mean:.3f}, "
        f"min={ppb_min:.3f}, "
        f"max={ppb_max:.3f}, "
        f"median={ppb_med:.3f}, "
        f"n={len(pp_bvals)}\n"
    )

    lf.write(
        "PP-PLDDT_CA:         "
        f"mean={df_mean:.3f}, "
        f"min={df_min:.3f}, "
        f"max={df_max:.3f}, "
        f"median={df_med:.3f}, "
        f"n={len(diff_vals)}\n"
    )

    lf.write(
        "PPxPLDDT_CA:         "
        f"mean={pr_mean:.3f}, "
        f"min={pr_min:.3f}, "
        f"max={pr_max:.3f}, "
        f"median={pr_med:.3f}, "
        f"n={len(prod_vals)}\n"
    )

print(
    f"[INFO] Site-wise PP statistics use the 0-100 scale "
    f"(n={len(pp_vals)} alignment columns)."
)
print(
    f"[INFO] Structure-based statistics use CA atoms "
    f"(n={len(pl_vals)} retained residues)."
)
print(f"[INFO] Log saved: {log_file}")
PY

echo "[INFO] --- Done ---"
echo "[INFO] Selected AlphaFold JSON: ${BEST_JSON}"
echo "[INFO] Selected AlphaFold CIF : ${BEST_CIF}"
echo "[INFO] Aligned with-gap FASTA : ${NODE_ALN_FASTA}"
echo "[INFO] Site-wise PP values    : ${TMP_PP_FILE}"
echo "[INFO] PDB with pLDDT         : ${PLDDT_PDB}"
echo "[INFO] PDB with PP            : ${PP_PDB}"
echo "[INFO] PDB with PP-pLDDT      : ${DIFF_PDB}"
echo "[INFO] PDB with PPxPLDDT      : ${PROD_PDB}"
echo "[INFO] Stats log              : ${STATS_LOG}"
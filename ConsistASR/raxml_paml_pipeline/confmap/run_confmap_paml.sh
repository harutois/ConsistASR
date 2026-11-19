#!/usr/bin/env bash
# ============================================================
# run_confmap_paml.sh
#
# Map site-wise posterior probabilities (PP) from a PAML rst
# file onto an AlphaFold model as PDB B-factors, and generate
# additional PDB files containing PP×pLDDT and scaled
# PP−pLDDT metrics.
#
# Workflow:
#   1. Detect the AlphaFold JSON/CIF pairs in the current
#      AlphaFold output directory and select the highest-ranked
#      valid model.
#   2. Convert the selected CIF to PDB, retaining AlphaFold
#      pLDDT values in the B-factor field.
#   3. Extract per-alignment-column PP values for the requested
#      PAML ancestral node from the rst marginal reconstruction
#      table using extract_pp_from_paml_rst.py.
#   4. Extract the corresponding gap-containing ancestral
#      sequence from the supplied with-gap FASTA.
#   5. Validate alignment-column count and exact sequence
#      correspondence between the ungapped ancestral sequence
#      and the AlphaFold structure.
#   6. Map PP onto structure residues using
#      map_confidence_to_bfactor.py.
#   7. Generate PDB files containing scaled PP−pLDDT and
#      PP×pLDDT values, matching residues explicitly rather
#      than relying on atom-list order.
#   8. Compute site-wise and CA-based confidence statistics.
#
# Required helper scripts (same directory as this wrapper):
#   - extract_pp_from_paml_rst.py
#   - map_confidence_to_bfactor.py
#
# Usage:
#   bash run_confmap_paml.sh \
#       --rst     <PAML_RST> \
#       --node    <PAML_NODE_INT> \
#       --withgap <WITHGAP_FASTA> \
#       [--outdir confmap]
#
# Example:
#   bash run_confmap_paml.sh \
#       --rst     rst \
#       --node    28 \
#       --withgap ASR_toy_7tm_rhodopsin_indel_withgap.fasta \
#       --outdir  confmap
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONSISTASR_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
COMMON_DIR="${CONSISTASR_ROOT}/common"

EXTRACT_PP_SCRIPT="${SCRIPT_DIR}/extract_pp_from_paml_rst.py"

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
RST_FILE=""
PAML_NODE=""
WITHGAP_FASTA=""
OUTDIR="confmap"

usage() {
    cat <<EOF
Usage:
  $0 --rst <PAML_RST> --node <PAML_NODE_INT> --withgap <WITHGAP_FASTA> [--outdir confmap]

Required:
  --rst PATH        PAML rst file
  --node INT        PAML ancestral-node ID (e.g. 28)
  --withgap PATH    Gap-containing ancestral FASTA preserving the
                    alignment coordinates represented in the PAML rst file

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
        --rst)
            require_value "$1" "${2-}"
            RST_FILE="$2"
            shift 2
            ;;
        --node)
            require_value "$1" "${2-}"
            PAML_NODE="$2"
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

if [[ -z "${RST_FILE}" || -z "${PAML_NODE}" || -z "${WITHGAP_FASTA}" ]]; then
    echo "[ERROR] --rst, --node, and --withgap are required." >&2
    usage >&2
    exit 1
fi

# -------------------------
# Basic checks
# -------------------------
if [[ ! -f "${RST_FILE}" ]]; then
    echo "[ERROR] RST_FILE not found: ${RST_FILE}" >&2
    exit 1
fi

if [[ ! -f "${WITHGAP_FASTA}" ]]; then
    echo "[ERROR] WITHGAP_FASTA not found: ${WITHGAP_FASTA}" >&2
    exit 1
fi

if [[ ! -f "${EXTRACT_PP_SCRIPT}" ]]; then
    echo "[ERROR] Required helper script not found: ${EXTRACT_PP_SCRIPT}" >&2
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

if [[ ! "${PAML_NODE}" =~ ^[0-9]+$ ]] || [[ "${PAML_NODE}" -le 0 ]]; then
    echo "[ERROR] --node must be a positive integer PAML node ID: ${PAML_NODE}" >&2
    exit 1
fi

NODE_NUM="${PAML_NODE}"
NODE_LABEL="Node${NODE_NUM}"

# ------------------------------------------------------------
# Create output directory
# (Assumes execution inside one AlphaFold Server job directory.)
# ------------------------------------------------------------
mkdir -p "${OUTDIR}"
OUTDIR_ABS="$(cd "${OUTDIR}" && pwd)"

echo "[INFO] Output directory: ${OUTDIR_ABS}"
echo "[INFO] PAML node: ${PAML_NODE} (label=${NODE_LABEL})"

PLDDT_PDB="${OUTDIR}/${NODE_LABEL}_plddt_bfactor.pdb"
TMP_PP_FILE="${OUTDIR}/${NODE_LABEL}_pp_values.txt"
NODE_ALN_FASTA="${OUTDIR}/${NODE_LABEL}_aligned_withgap.fasta"
PP_PDB="${OUTDIR}/${NODE_LABEL}_pp_bfactor.pdb"
DIFF_PDB="${OUTDIR}/${NODE_LABEL}_ppminusplddt_bfactor.pdb"
PROD_PDB="${OUTDIR}/${NODE_LABEL}_ppxplddt_bfactor.pdb"
STATS_LOG="${OUTDIR}/${NODE_LABEL}_conf_CA_stats.log"

# ------------------------------------------------------------
# Step 1) Detect AlphaFold JSON/CIF pairs
#
# The user-defined AlphaFold Server job name is NOT used to
# identify the PAML node. The current directory is expected to
# contain output from one AlphaFold job, and exact sequence
# correspondence is validated downstream.
# ------------------------------------------------------------
echo "[INFO] --- Step 1: Searching AlphaFold model outputs ---"

BEST_INFO=$("${PYTHON_BIN}" - <<'PY'
import glob
import json
import math
import os
import re
import sys

candidates = sorted(
    glob.glob("*_summary_confidences_*.json")
)

if not candidates:
    sys.exit(
        "[ERROR] No *_summary_confidences_*.json found "
        "in the current directory."
    )

parsed = []

for path in candidates:
    m = re.match(
        r"^(.*)_summary_confidences_(\d+)\.json$",
        os.path.basename(path),
    )

    if not m:
        print(
            f"[WARN] Skipping unrecognized summary filename: "
            f"{path}",
            file=sys.stderr,
        )
        continue

    prefix = m.group(1)
    index = int(m.group(2))

    parsed.append(
        (path, prefix, index)
    )

if not parsed:
    sys.exit(
        "[ERROR] No summary-confidence filename "
        "could be parsed."
    )

prefixes = sorted(
    {prefix for _, prefix, _ in parsed}
)

if len(prefixes) != 1:
    msg = [
        "[ERROR] Multiple AlphaFold output prefixes were found "
        "in the current directory.",
        "The confidence-mapping wrapper expects outputs from "
        "one AlphaFold job per working directory.",
        "Detected prefixes:",
    ]

    msg.extend(
        f"  - {prefix}"
        for prefix in prefixes
    )

    sys.exit("\n".join(msg))

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
            f"[WARN] Skipping unreadable JSON "
            f"{json_path}: {exc}",
            file=sys.stderr,
        )
        continue

    score_value = data.get("ranking_score")

    if score_value is None:
        score_value = data.get("plddt")

    if score_value is None:
        print(
            f"[WARN] Skipping JSON without ranking_score "
            f"or plddt: {json_path}",
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
        "[ERROR] No valid AlphaFold summary-confidence "
        "JSON/CIF pair was found."
    )

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
    echo "[ERROR] Failed to select an AlphaFold model." >&2
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

# ------------------------------------------------------------
# Step 2) CIF -> PDB conversion
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
# Step 3) Extract PP values from PAML rst
# ------------------------------------------------------------
echo "[INFO] --- Step 3: Extracting PP from PAML RST for node ${PAML_NODE} ---"

"${PYTHON_BIN}" "${EXTRACT_PP_SCRIPT}" \
    --rst "${RST_FILE}" \
    --paml_node "${PAML_NODE}" \
    --withgap_fasta "${WITHGAP_FASTA}" \
    --out "${TMP_PP_FILE}"

if [[ ! -s "${TMP_PP_FILE}" ]]; then
    echo "[ERROR] PP file was not generated: ${TMP_PP_FILE}" >&2
    exit 1
fi

# ------------------------------------------------------------
# Step 4) Extract exact ancestral sequence from with-gap FASTA
# ------------------------------------------------------------
echo "[INFO] --- Step 4: Extracting aligned sequence for node ${PAML_NODE} ---"

"${PYTHON_BIN}" - \
    "${WITHGAP_FASTA}" \
    "${PAML_NODE}" \
    "${NODE_ALN_FASTA}" <<'PY'
import sys

fasta, node, out_file = sys.argv[1:]

matches = []
header = None
seq_parts = []


def store_current():
    if header == node:
        matches.append(
            "".join(seq_parts)
        )


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
                    f"[ERROR] Sequence data found before "
                    f"the first FASTA header at line {lineno} "
                    f"in {fasta}."
                )

            seq_parts.append(
                "".join(line.split())
            )

store_current()

if len(matches) == 0:
    sys.exit(
        f"[ERROR] Exact FASTA header >{node} "
        f"not found in {fasta}."
    )

if len(matches) > 1:
    sys.exit(
        f"[ERROR] Duplicate exact FASTA header >{node} "
        f"found in {fasta}."
    )

sequence = matches[0]

if not sequence:
    sys.exit(
        f"[ERROR] Empty sequence found for >{node} "
        f"in {fasta}."
    )

with open(out_file, "w") as out:
    out.write(f">{node}\n")
    out.write(sequence + "\n")

print(
    f"[INFO] Extracted aligned sequence for node {node}: "
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
    "${NODE_LABEL}" <<'PY'
import sys

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

fasta_file, pp_file, pdb_file, node_label = sys.argv[1:]

# Read aligned ancestral sequence
seq_parts = []

with open(fasta_file) as fh:
    for raw in fh:
        line = raw.strip()

        if not line or line.startswith(">"):
            continue

        seq_parts.append(line)

aligned = "".join(seq_parts).upper()

if not aligned:
    sys.exit(
        "[ERROR] Extracted ancestral sequence is empty."
    )

ungapped = aligned.replace("-", "")

# Read PP values
pp_vals = []

with open(pp_file) as fh:
    for lineno, raw in enumerate(fh, start=1):
        value = raw.strip()

        if not value:
            continue

        try:
            pp_vals.append(float(value))
        except ValueError:
            sys.exit(
                f"[ERROR] Invalid PP value at line {lineno} "
                f"in {pp_file}: {value!r}"
            )

if len(pp_vals) != len(aligned):
    sys.exit(
        f"[ERROR] Alignment-column mismatch for {node_label}: "
        f"{len(pp_vals)} PP sites versus "
        f"{len(aligned)} columns in the with-gap "
        f"ancestral sequence."
    )

# Extract protein sequence from converted AlphaFold PDB
parser = PDBParser(QUIET=True)

structure = parser.get_structure(
    "AF",
    pdb_file,
)

models = list(structure)

if len(models) != 1:
    sys.exit(
        f"[ERROR] Expected a single AlphaFold model in "
        f"{pdb_file}, but found {len(models)} models."
    )

protein_chains = []

for chain in models[0]:
    residues = []

    for residue in chain:
        if not is_aa(
            residue,
            standard=False,
        ):
            continue

        if "CA" not in residue:
            sys.exit(
                f"[ERROR] Protein residue "
                f"{residue.get_resname()} {residue.id[1]} "
                f"in chain {chain.id!r} has no CA atom."
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
        f"but found {len(protein_chains)} protein chains "
        f"({desc})."
    )

chain_id, structure_seq = protein_chains[0]

if len(ungapped) != len(structure_seq):
    sys.exit(
        f"[ERROR] Ungapped ancestral length != "
        f"structure sequence length for {node_label}: "
        f"ancestral={len(ungapped)}, "
        f"structure={len(structure_seq)}."
    )

if ungapped != structure_seq:
    mismatch = None

    for i, (anc_aa, struct_aa) in enumerate(
        zip(ungapped, structure_seq),
        start=1,
    ):
        if anc_aa != struct_aa:
            mismatch = (
                i,
                anc_aa,
                struct_aa,
            )
            break

    if mismatch is None:
        detail = "sequence mismatch detected"
    else:
        pos, anc_aa, struct_aa = mismatch

        detail = (
            f"first mismatch at residue {pos}: "
            f"ancestral={anc_aa}, "
            f"structure={struct_aa}"
        )

    sys.exit(
        f"[ERROR] AlphaFold structure sequence does not "
        f"exactly match the ungapped ancestral sequence "
        f"for {node_label}; {detail}."
    )

print(
    f"[INFO] PP/alignment columns: "
    f"{len(pp_vals)} / {len(aligned)}"
)

print(
    f"[INFO] Ungapped ancestral sequence: "
    f"{len(ungapped)} residues"
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
# Step 6) Embed PP into B-factors
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
# Step 7) Generate PP-pLDDT and PPxpLDDT PDBs
#
# Residues are matched explicitly by model / chain / residue
# identifier / residue name. Atom-list order is not used.
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

pl_structure = parser.get_structure(
    "PL",
    plddt_pdb,
)

pp_structure = parser.get_structure(
    "PP",
    pp_pdb,
)


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

                key = residue_key(
                    model,
                    chain,
                    residue,
                )

                if key in values:
                    raise RuntimeError(
                        f"Duplicate residue key in {label}: {key}"
                    )

                values[key] = (
                    residue["CA"].get_bfactor()
                )

    if not values:
        raise RuntimeError(
            f"No CA-containing residues found in {label}."
        )

    return values


try:
    pl_map = ca_bfactor_map(
        pl_structure,
        "pLDDT PDB",
    )

    pp_map = ca_bfactor_map(
        pp_structure,
        "PP PDB",
    )

except RuntimeError as exc:
    sys.exit(f"[ERROR] {exc}")

pl_keys = set(pl_map)
pp_keys = set(pp_map)

if pl_keys != pp_keys:
    only_pl = sorted(
        pl_keys - pp_keys,
        key=str,
    )

    only_pp = sorted(
        pp_keys - pl_keys,
        key=str,
    )

    msg = [
        "[ERROR] Residue identifiers differ between "
        "the pLDDT and PP PDB files."
    ]

    if only_pl:
        msg.append(
            "Present only in pLDDT PDB: "
            + ", ".join(
                map(str, only_pl[:5])
            )
        )

    if only_pp:
        msg.append(
            "Present only in PP PDB: "
            + ", ".join(
                map(str, only_pp[:5])
            )
        )

    sys.exit("\n".join(msg))

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

    # PP - pLDDT: [-100, +100]
    diff = pp - pl

    # Scaled to [0, 100]
    scaled_diff = (
        diff + 100.0
    ) / 2.0

    # PP x pLDDT / 100
    prod = (
        pp * pl
    ) / 100.0

    diff_values[key] = scaled_diff
    prod_values[key] = prod


def apply_values(structure, values):
    seen = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" not in residue:
                    continue

                key = residue_key(
                    model,
                    chain,
                    residue,
                )

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


diff_structure = deepcopy(
    pp_structure
)

prod_structure = deepcopy(
    pp_structure
)

try:
    apply_values(
        diff_structure,
        diff_values,
    )

    apply_values(
        prod_structure,
        prod_values,
    )

except RuntimeError as exc:
    sys.exit(f"[ERROR] {exc}")

io = PDBIO()

io.set_structure(
    diff_structure
)
io.save(
    diff_pdb
)

io.set_structure(
    prod_structure
)
io.save(
    prod_pdb
)

print(
    f"[INFO] Matched {len(pp_keys)} residues "
    "between PP and pLDDT structures."
)

print(
    f"[INFO] Wrote diff PDB "
    f"(PP - pLDDT scaled): {diff_pdb}"
)

print(
    f"[INFO] Wrote prod PDB "
    f"(PP x pLDDT):       {prod_pdb}"
)
PY

# ------------------------------------------------------------
# Step 8) Statistics
# ------------------------------------------------------------
echo "[INFO] --- Step 8: Computing confidence statistics ---"

"${PYTHON_BIN}" - \
    "${TMP_PP_FILE}" \
    "${PLDDT_PDB}" \
    "${PP_PDB}" \
    "${DIFF_PDB}" \
    "${PROD_PDB}" \
    "${STATS_LOG}" \
    "${NODE_LABEL}" <<'PY'
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
    node_label,
) = sys.argv[1:]


def stats(arr, label):
    if not arr:
        sys.exit(
            f"[ERROR] No values available for statistics: "
            f"{label}"
        )

    return (
        statistics.mean(arr),
        min(arr),
        max(arr),
        statistics.median(arr),
    )


# -------------------------
# Site-wise PP from rst
# -------------------------
pp_vals_raw = []

with open(pp_file) as fh:
    for lineno, raw in enumerate(fh, start=1):
        value = raw.strip()

        if not value:
            continue

        try:
            pp_vals_raw.append(
                float(value)
            )
        except ValueError:
            sys.exit(
                f"[ERROR] Invalid PP value at line {lineno} "
                f"in {pp_file}: {value!r}"
            )

if not pp_vals_raw:
    sys.exit(
        "[ERROR] No site-wise PP values were loaded."
    )

if min(pp_vals_raw) < 0:
    sys.exit(
        "[ERROR] Negative PP value detected."
    )

# PAML rst PP is normally on the 0-1 scale.
# Retain compatibility with an already 0-100 input.
if max(pp_vals_raw) <= 1.0:
    pp_vals = [
        value * 100.0
        for value in pp_vals_raw
    ]

elif max(pp_vals_raw) <= 100.0:
    pp_vals = pp_vals_raw

else:
    sys.exit(
        "[ERROR] Site-wise PP values exceed "
        "the expected 0-100 range."
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
    structure = parser.get_structure(
        path,
        path,
    )

    return [
        atom.get_bfactor()
        for atom in structure.get_atoms()
        if atom.get_name() == "CA"
    ]


pl_vals = bfacts_ca(
    plddt_pdb
)

pp_bvals = bfacts_ca(
    pp_pdb
)

diff_vals = bfacts_ca(
    diff_pdb
)

prod_vals = bfacts_ca(
    prod_pdb
)

n_ca = {
    "pLDDT": len(pl_vals),
    "PP": len(pp_bvals),
    "PP-pLDDT": len(diff_vals),
    "PPxpLDDT": len(prod_vals),
}

if len(set(n_ca.values())) != 1:
    sys.exit(
        "[ERROR] CA counts differ among confidence-mapped "
        "PDB files: "
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
        f"[INFO] Node {node_label} confidence statistics\n"
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
    f"[INFO] Site-wise PP statistics use the "
    f"0-100 scale (n={len(pp_vals)} alignment columns)."
)

print(
    f"[INFO] Structure-based statistics use CA atoms "
    f"(n={len(pl_vals)} retained residues)."
)

print(
    f"[INFO] Log saved: {log_file}"
)
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
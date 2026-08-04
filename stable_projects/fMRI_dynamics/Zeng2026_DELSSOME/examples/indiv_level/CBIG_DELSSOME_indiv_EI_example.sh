#!/usr/bin/env bash
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# Co-authored by Fang Tian

# Individual-level DELSSOME demo: FC + FCD-CDF CSVs in, E/I ratio CSV out.
#
# This is the whole user-facing workflow in one script:
#   1. CBIG_DELSSOME_indiv_estimate_EI.py  -- once per subject
#   2. CBIG_DELSSOME_indiv_collect_EI.py   -- assemble the GAMLSS input table
#   3. the checker                         -- diff against reference_output/
#
# The bundled data is SYNTHETIC (see ../../README.md, "Example data"), so the
# E/I values it produces are not a scientific result -- they only demonstrate
# that the pipeline runs and reproduces bit-for-bit.
#
# Usage:
#     bash examples/indiv_level/CBIG_DELSSOME_indiv_EI_example.sh
#
# Overridable: SUB_IDS, NUM_EPOCHS, SEED, DELSSOME_CKPT, SKIP_CHECK

set -euo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO="$( cd "${HERE}/../.." && pwd )"
INDIV="${REPO}/indiv_level"

cd "${INDIV}"
export PYTHONPATH="${INDIV}:${PYTHONPATH:-}"
# Keep the dataset cache out of any read-only input mount.
export DELSSOME_TMP_DIR="${HERE}/output/tmp"

INPUT="${HERE}/input"
OUTPUT="${HERE}/output"
CONFIG="${HERE}/config/example_indiv_pFIC.ini"

SUB_IDS="${SUB_IDS:-sub-000 sub-001}"
# Tiny epoch count so the demo finishes in minutes rather than hours.
NUM_EPOCHS="${NUM_EPOCHS:-3}"
SEED="${SEED:-123456}"
# The paper's individual-level predictor, if it has been downloaded. Optional:
# without it the demo runs pure Euler, which needs no predictor at all.
DELSSOME_CKPT="${DELSSOME_CKPT:-${REPO}/pretrained_models/indiv_level/pFIC_DK68_predictor.ckpt}"

# TR and scan length of the (synthetic) demo scans. There is deliberately no
# default for these in estimate_EI: silently inheriting HCP-YA's 0.72 s would
# produce wrong numbers for anyone else's data.
TR="${TR:-0.72}"
# 60 s rather than the real 864 s, purely for speed.
SCAN_LENGTH="${SCAN_LENGTH:-60}"

mkdir -p "${OUTPUT}"

# --- 1. one E/I estimate per subject (Euler) -------------------------------- #
for sub_id in ${SUB_IDS}; do
    echo ""
    echo "======== Euler CMA-ES :: ${sub_id} ========"
    python -u -m scripts.CBIG_DELSSOME_indiv_estimate_EI \
        --fc          "${INPUT}/${sub_id}_FC.csv" \
        --fcd-cdf     "${INPUT}/${sub_id}_FCD_CDF.csv" \
        --sc          "${INPUT}/SC_group_DK68.csv" \
        --myelin      "${INPUT}/myelin_DK68.csv" \
        --gradient    "${INPUT}/rsfc_gradient_DK68.csv" \
        --tr          "${TR}" \
        --scan-length "${SCAN_LENGTH}" \
        --config      "${CONFIG}" \
        --num-epochs  "${NUM_EPOCHS}" \
        --seed        "${SEED}" \
        --subject-id  "${sub_id}" \
        --overwrite \
        --out-dir     "${OUTPUT}/euler/${sub_id}"
done

# --- 2. assemble the GAMLSS input table ------------------------------------ #
echo ""
echo "======== collect E/I ratios into a GAMLSS table ========"
python -u -m scripts.CBIG_DELSSOME_indiv_collect_EI \
    --in-dir "${OUTPUT}/euler" \
    --demogr "${INPUT}/demogr.csv" \
    --out    "${OUTPUT}/euler/gamlss_input_from_demo.csv"

# --- 3. the same thing with DELSSOME acceleration -------------------------- #
# Identical command plus --delssome-ckpt. Only the SEARCH stage changes engine;
# validation and test stay on Euler, following the paper.
if [ -f "${DELSSOME_CKPT}" ]; then
    first_sub="${SUB_IDS%% *}"
    echo ""
    echo "======== DELSSOME-accelerated CMA-ES :: ${first_sub} ========"
    python -u -m scripts.CBIG_DELSSOME_indiv_estimate_EI \
        --fc            "${INPUT}/${first_sub}_FC.csv" \
        --fcd-cdf       "${INPUT}/${first_sub}_FCD_CDF.csv" \
        --sc            "${INPUT}/SC_group_DK68.csv" \
        --myelin        "${INPUT}/myelin_DK68.csv" \
        --gradient      "${INPUT}/rsfc_gradient_DK68.csv" \
        --tr            "${TR}" \
        --scan-length   "${SCAN_LENGTH}" \
        --config        "${CONFIG}" \
        --num-epochs    "${NUM_EPOCHS}" \
        --seed          "${SEED}" \
        --subject-id    "${first_sub}" \
        --delssome-ckpt "${DELSSOME_CKPT}" \
        --overwrite \
        --out-dir       "${OUTPUT}/delssome/${first_sub}"
else
    echo ""
    echo "[skip] DELSSOME-accelerated run: no predictor at"
    echo "         ${DELSSOME_CKPT}"
    echo "       The Euler runs above already produced E/I ratios -- the"
    echo "       predictor only makes the search stage cheaper to evaluate."
    echo "       See ../../README.md"
    echo "       (\"Data\") for where to obtain it."
fi

# --- 4. regression check --------------------------------------------------- #
if [ "${SKIP_CHECK:-0}" != "1" ]; then
    echo ""
    echo "======== checking against reference_output/ ========"
    python -u "${HERE}/CBIG_DELSSOME_indiv_check_EI_example_results.py"
fi

echo ""
echo "Done. Results:"
echo "  ${OUTPUT}/euler/<sub_id>/EI_ratio.csv          per-region E/I ratio"
echo "  ${OUTPUT}/euler/<sub_id>/best_parameters.csv    fitted wEE / wEI / sigma / G"
echo "  ${OUTPUT}/euler/<sub_id>/fit_quality.csv        FC+FCD cost per stage"
echo "  ${OUTPUT}/euler/gamlss_input_from_demo.csv      input for the GAMLSS demo"

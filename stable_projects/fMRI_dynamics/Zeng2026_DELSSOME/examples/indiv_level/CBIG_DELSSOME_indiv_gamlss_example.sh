#!/usr/bin/env bash
# Written by Fang Tian and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# Co-authored by Tianchu Zeng

# GAMLSS demo: a table of per-subject E/I ratios in, a normative lifespan
# trajectory out, then the fitted model applied to a second cohort to show the
# site-harmonisation workflow.
#
# Note this uses a bundled 600-subject SYNTHETIC table, NOT the output of
# CBIG_DELSSOME_indiv_EI_example.sh. That demo fits two subjects; GAMLSS needs
# several hundred before a normative curve means anything. Run
# CBIG_DELSSOME_indiv_collect_EI.py on your own cohort to build the real table
# -- input/gamlss_input.csv shows the exact format it produces.
#
# Usage:
#     bash examples/indiv_level/CBIG_DELSSOME_indiv_gamlss_example.sh
#
# Overridable: CONFIG_NAME (any JSON under indiv_level/GAMLSS/configs/), SKIP_CHECK

set -euo pipefail
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO="$( cd "${HERE}/../.." && pwd )"
GAMLSS_DIR="${REPO}/indiv_level/GAMLSS"

INPUT="${HERE}/input"
OUTPUT="${HERE}/output/gamlss"
CONFIG_NAME="${CONFIG_NAME:-SHASHo2_bs_df4.json}"
CONFIG_STEM="${CONFIG_NAME%.json}"

if ! command -v Rscript >/dev/null 2>&1; then
    echo "[skip] Rscript not available. Install R + gamlss (see"
    echo "       ../../indiv_level/GAMLSS/README.md) and re-run."
    exit 0
fi

for csv in gamlss_input.csv gamlss_new_cohort.csv; do
    if [ ! -f "${INPUT}/${csv}" ]; then
        echo "[ERROR] Missing bundled input ${INPUT}/${csv}"
        exit 1
    fi
done

mkdir -p "${OUTPUT}"

# Both scripts take absolute paths and locate their own helpers, so there is no
# need to cd into the source tree or stage anything into it.
echo "1/2 Fitting GAMLSS with configs/${CONFIG_NAME} ..."
Rscript "${GAMLSS_DIR}/CBIG_fit_gamlss_model.R" \
    "${GAMLSS_DIR}/configs/${CONFIG_NAME}" \
    --input_data "${INPUT}/gamlss_input.csv" \
    --output_dir "${OUTPUT}/${CONFIG_STEM}"

echo "2/2 Applying the fitted model to the bundled new-cohort table ..."
Rscript "${GAMLSS_DIR}/CBIG_apply_gamlss_model.R" \
    "${OUTPUT}/${CONFIG_STEM}/fitted_model.rds" \
    "${INPUT}/gamlss_new_cohort.csv" \
    "${OUTPUT}/apply_demo/${CONFIG_STEM}_outputs.csv"

if [ "${SKIP_CHECK:-0}" != "1" ]; then
    echo ""
    echo "======== checking against reference_output/ ========"
    CONFIG_NAME="${CONFIG_NAME}" \
        python -u "${HERE}/CBIG_DELSSOME_indiv_check_gamlss_example_results.py"
fi

echo ""
echo "Done. Results:"
echo "  ${OUTPUT}/${CONFIG_STEM}/            fitted model, centiles, summary"
echo "  ${OUTPUT}/apply_demo/                harmonised values for the new cohort"

#!/bin/sh
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Generate a stand-alone copy of the Zeng2026_DELSSOME project (plus shared
# dependencies) by copying the whole CBIG repository and removing the stable
# projects that Zeng2026_DELSSOME does not depend on.

# CBIG_CODE_DIR is set by CBIG_DELSSOME_tested_config.sh; source that first.
if [ -z "$CBIG_CODE_DIR" ]; then
    echo "[error] CBIG_CODE_DIR is not set. Source CBIG_DELSSOME_tested_config.sh first." >&2
    exit 1
fi

cd "$(dirname "$CBIG_CODE_DIR")"
rsync -a --exclude .git "$CBIG_CODE_DIR"/* Standalone_Zeng2026_DELSSOME

# remove unused stable projects
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/brain_parcellation
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/disorder_subtypes
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/harmonization
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/meta-analysis
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/predict_phenotypes
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/preprocessing
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/registration

# remove sibling fMRI_dynamics projects
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/fMRI_dynamics/Kong2021_pMFM
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/fMRI_dynamics/Liegeois2017_Surrogates
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/fMRI_dynamics/Wang2018_MFMem
rm -r Standalone_Zeng2026_DELSSOME/stable_projects/fMRI_dynamics/Zhang2024_pFIC

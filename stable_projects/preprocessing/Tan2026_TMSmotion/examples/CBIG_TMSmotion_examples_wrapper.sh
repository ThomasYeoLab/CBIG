#!/bin/bash
# Written by Trevor Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#
# This examples wrapper script runs and checks simulations of rs-fMRI.
#
# Example:
# ${TAN2026_TMSMOTION_DIR}/examples/CBIG_TMSmotion_examples_wrapper.sh;

# -------------------------------------------------------------------
# Check that TAN2026_TMSMOTION_DIR is set
# -------------------------------------------------------------------
if [ -z "$TAN2026_TMSMOTION_DIR" ]; then
    echo "ERROR: TAN2026_TMSMOTION_DIR environment variable is not set."
    echo "Please set it before running this script, e.g.:"
    echo "  export TAN2026_TMSMOTION_DIR=/path/to/Tan2026_TMSmotion"
    exit 1
fi

# -------------------------------------------------------------------
# Define script paths
# -------------------------------------------------------------------
EXAMPLES_SCRIPT="$TAN2026_TMSMOTION_DIR/examples/scripts/CBIG_TMSmotion_examples.py"
CHECK_SCRIPT="$TAN2026_TMSMOTION_DIR/examples/scripts/CBIG_TMSmotion_check_examples_result.py"
source CBIG_init_conda; conda activate CBIG_TMSmotion;

# -------------------------------------------------------------------
# Run examples simulation
# -------------------------------------------------------------------
echo "Running example simulations..."
python "$EXAMPLES_SCRIPT" "examples"

# -------------------------------------------------------------------
# Check example results
# -------------------------------------------------------------------
echo "Checking example results..."
python "$CHECK_SCRIPT"
conda deactivate;

echo "All example scripts completed successfully."
#!/usr/bin/env bash
# Written by Trevor Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#
# All-in-one driver for the Tan2026_TMSmotion replication pipeline, for 1
# pipeline (pl_name) at a time. Runs in sequence:
#   init_project -> init_pipeline -> select_subID_window -> select_subID_sim
#   -> (TMS_Kong2025 gt & ngt in parallel) -> compare_Kong2025
#   -> check_reference_output
#
# gt & ngt are run in parallel since they're independent and compare_Kong2025
# needs both to finish anyway. CBIG_TMSmotion_check_reference_output.py runs directly on
# this terminal (it's a lightweight comparison, not a cluster job).
#
# The terminal only shows which step is currently running and whether it
# passed/failed -- not each step's full output. Full output for every step
# always goes to its own per-step log file under
# replication/logs/<pl_name>_<mode>_<timestamp>/. The step status lines
# shown on the terminal are also mirrored into one consolidated log file:
#   replication/<pl_name>/output/replication/replication_<pl_name>_<timestamp>.log
# On any failure, the terminal prints which step failed and reminds you to
# check the log folder.
#
# To run 2 pipelines (e.g. pl1 & pl2) at once, invoke this script once per
# pipeline in separate terminals -- job names already include pl_name, so
# concurrent invocations for different pipelines don't interfere with each
# other's `qstat` polling, and each pipeline's outputs live under its own
# replication/<pl_name>/output/ directory.
#
# This must be run on headnode (it submits jobs via qsub/CBIG_pbsubmit).
#
# Usage:
#   CBIG_TMSmotion_replication_wrapper.sh <full|lite> <pl_name>
#
# Example:
#   CBIG_TMSmotion_replication_wrapper.sh lite pl1

set -uo pipefail

#############################################################################
# ARGUMENT PARSING
#############################################################################
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <full|lite> <pl_name>" >&2
    exit 1
fi

MODE="$1"
PL_NAME="$2"

if [[ "$MODE" != "full" && "$MODE" != "lite" ]]; then
    echo "ERROR: mode must be 'full' or 'lite', got '$MODE'" >&2
    exit 1
fi

#############################################################################
# ENVIRONMENT CHECKS
#############################################################################
if [[ -z "${TAN2026_TMSMOTION_DIR:-}" ]]; then
    echo "ERROR: TAN2026_TMSMOTION_DIR is not set." >&2
    exit 1
fi

if [[ -z "${CBIG_REPDATA_DIR:-}" ]]; then
    echo "ERROR: CBIG_REPDATA_DIR is not set." >&2
    exit 1
fi

CURRNODE=$(hostname)
if [[ "$CURRNODE" != "headnode" ]]; then
    echo "ERROR: this wrapper submits jobs and must be run on headnode (currently on '$CURRNODE')." >&2
    exit 1
fi

REP_DIR="$TAN2026_TMSMOTION_DIR/replication"
INI_FILE="$REP_DIR/init/input/${PL_NAME}.ini"
PL_CONFIG="$REP_DIR/${PL_NAME}/config/pipeline.json"

if [[ ! -f "$INI_FILE" ]]; then
    echo "ERROR: ini file not found: $INI_FILE" >&2
    exit 1
fi

#############################################################################
# LOGGING
#############################################################################
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="$REP_DIR/logs/${PL_NAME}_${MODE}_${TIMESTAMP}"
mkdir -p "$LOG_DIR"

# Consolidated whole-run log: mirrors whatever this script prints to the
# terminal (step status lines, not each step's full output -- that goes
# only to its own per-step log file under LOG_DIR, see run_step below).
REPLICATION_LOG_DIR="$REP_DIR/${PL_NAME}/output/replication"
mkdir -p "$REPLICATION_LOG_DIR"
REPLICATION_LOG="$REPLICATION_LOG_DIR/replication_${PL_NAME}_${TIMESTAMP}.log"
exec > >(tee -a "$REPLICATION_LOG") 2>&1

echo "Mode: $MODE | Pipeline: $PL_NAME"
echo "Logs: $LOG_DIR"
echo "Replication log: $REPLICATION_LOG"

# Some steps (e.g. `03_select_subID_sim` in full mode, TMS_Kong2025) can run
# for hours with no output in between, which can cause the headnode session
# to be logged out for inactivity. HEARTBEAT_INTERVAL_SEC controls how often
# a "still running" line is printed while a step's real command runs, purely
# to keep the terminal/SSH connection from looking idle.
HEARTBEAT_INTERVAL_SEC=60

# run_step <log_name> <cmd...>
# Prints a concise "running/done/failed" status to the terminal (and the
# consolidated log), plus a periodic heartbeat line while cmd is still
# running. cmd itself still runs in the foreground exactly as before (only
# a background heartbeat printer is added alongside it, killed as soon as
# cmd finishes) -- its full stdout+stderr goes ONLY to its own per-step log
# file, not to the terminal. Combined with `set -e` below, a failing step
# aborts the whole script.
run_step() {
    local log_name="$1"
    shift
    local log_file="$LOG_DIR/${log_name}.log"
    echo "Running step: $log_name ..."

    ( while true; do sleep "$HEARTBEAT_INTERVAL_SEC"; echo "  ... still running $log_name"; done ) &
    local heartbeat_pid=$!

    local had_e=0
    case "$-" in *e*) had_e=1 ;; esac
    set +e
    "$@" > "$log_file" 2>&1
    local status=$?
    [[ $had_e -eq 1 ]] && set -e

    kill "$heartbeat_pid" 2>/dev/null || true
    wait "$heartbeat_pid" 2>/dev/null || true

    if [[ $status -ne 0 ]]; then
        echo "ERROR: step '$log_name' failed (exit $status). See $log_file"
        return 1
    fi
    echo "Step '$log_name' completed."
    return 0
}

on_failure() {
    echo "Replication unsuccessful for $PL_NAME ($MODE)."
    echo "Please check log folder under: $LOG_DIR"
}
trap on_failure ERR

set -e

#############################################################################
# STEP: INIT PROJECT
#############################################################################
run_step "00_init_project" \
    python "$REP_DIR/init/script/CBIG_TMSmotion_init_project.py"

#############################################################################
# STEP: INIT PIPELINE
#############################################################################
run_step "01_init_pipeline" \
    python "$REP_DIR/init/script/CBIG_TMSmotion_init_pipeline.py" "$INI_FILE" "$MODE"

#############################################################################
# STEP: SELECT_SUBID_WINDOW
#############################################################################
run_step "02_select_subID_window" \
    python "$REP_DIR/select_subID_window/script/CBIG_TMSmotion_select_subID_window_wrapper.py" "$PL_CONFIG"

#############################################################################
# STEP: SELECT_SUBID_SIM
#############################################################################
run_step "03_select_subID_sim" \
    python "$REP_DIR/select_subID_sim/script/CBIG_TMSmotion_select_subID_sim_wrapper.py" "$PL_CONFIG"

#############################################################################
# STEP: TMS_KONG2025 GT & NGT, IN PARALLEL
#############################################################################
echo "Running step: 04_TMS_Kong2025_gt & 04_TMS_Kong2025_ngt (in parallel) ..."
GT_LOG="$LOG_DIR/04_TMS_Kong2025_gt.log"
NGT_LOG="$LOG_DIR/04_TMS_Kong2025_ngt.log"

python "$REP_DIR/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_gt_wrapper.py" "$PL_CONFIG" \
    > "$GT_LOG" 2>&1 &
GT_PID=$!

python "$REP_DIR/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_ngt_wrapper.py" "$PL_CONFIG" \
    > "$NGT_LOG" 2>&1 &
NGT_PID=$!

( while true; do sleep "$HEARTBEAT_INTERVAL_SEC"; echo "  ... still running 04_TMS_Kong2025_gt/ngt"; done ) &
HEARTBEAT_PID=$!

GT_STATUS=0
NGT_STATUS=0
wait "$GT_PID" || GT_STATUS=$?
wait "$NGT_PID" || NGT_STATUS=$?

kill "$HEARTBEAT_PID" 2>/dev/null || true
wait "$HEARTBEAT_PID" 2>/dev/null || true

if [[ $GT_STATUS -ne 0 ]]; then
    echo "ERROR: step '04_TMS_Kong2025_gt' failed (exit $GT_STATUS). See $GT_LOG"
fi
if [[ $NGT_STATUS -ne 0 ]]; then
    echo "ERROR: step '04_TMS_Kong2025_ngt' failed (exit $NGT_STATUS). See $NGT_LOG"
fi
if [[ $GT_STATUS -ne 0 || $NGT_STATUS -ne 0 ]]; then
    on_failure
    exit 1
fi
echo "Step '04_TMS_Kong2025_gt' completed."
echo "Step '04_TMS_Kong2025_ngt' completed."

#############################################################################
# STEP: COMPARE_KONG2025 (needs both gt & ngt to have succeeded)
#############################################################################
run_step "05_compare_Kong2025" \
    python "$REP_DIR/compare_Kong2025/script/CBIG_TMSmotion_compare_Kong2025_wrapper.py" "$PL_CONFIG"

#############################################################################
# STEP: CHECK REFERENCE OUTPUT (runs directly, no job submission)
#############################################################################
run_step "06_check_reference_output" \
    python "$REP_DIR/CBIG_TMSmotion_check_reference_output.py" "$PL_NAME" "$MODE"

trap - ERR
echo "Replication successful for $PL_NAME ($MODE)."
echo "Logs: $LOG_DIR"
echo "Replication log: $REPLICATION_LOG"

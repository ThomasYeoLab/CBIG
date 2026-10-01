'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Examples:


python $TAN2026_TMSMOTION_DIR/replication/init/script/CBIG_TMSmotion_init_pipeline.py \
    $TAN2026_TMSMOTION_DIR/replication/init/input/pl1.ini lite

python $TAN2026_TMSMOTION_DIR/replication/init/script/CBIG_TMSmotion_init_pipeline.py \
    $TAN2026_TMSMOTION_DIR/replication/init/input/pl2.ini lite

"""
# ===========================================================================
# IMPORTS
# ===========================================================================
import os
import sys
import configparser

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
CBIG_REPDATA_DIR = os.environ.get('CBIG_REPDATA_DIR')
import CBIG_TMSmotion_util as util

CBIG_CODE_DIR = util.get_cbig_code_dir()

sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'preproc', 'struct'))


# ===========================================================================
# FUNCTIONS
# ===========================================================================
def generate_pl_config(TAN2026_TMSMOTION_DIR, mode, **kwargs):
    """Build the pipeline configuration dictionary and write it to a JSON file.

    Assembles paths, preprocessing parameters, and subID-selection/
    simulation parameters (sourced from `**kwargs`, which originates from
    the input `.ini` file) into a single nested config dictionary, then
    writes it to
    `<TAN2026_TMSMOTION_DIR>/replication/<pl_name>/config/pipeline.json`.
    Sourcing parameters from a single `.ini` file (rather than editing
    this script directly) keeps changes to pipeline arguments centralized
    and the code clean.

    Args:
        TAN2026_TMSMOTION_DIR (str): Root directory of the
            Tan2026_TMSmotion project.
        mode (str): Either "full" or "lite". Selects the input data scale
            (input_data/<pl_name>/<mode>/...). Raw simulation generation
            (num_sim_per_bin etc.) always runs at full scale regardless
            of mode; lite-mode restriction to `lite_num_sim_ses`
            simulated sessions happens as a post-hoc filter downstream,
            in generate_select_subID_summary().
        **kwargs: Pipeline parameters read from the input `.ini` file
            (e.g. `FD_metric`, `DV_metric`, `window_time_sec`, etc.); see
            `read_pl_ini` for the full set of expected keys. Must include
            `pl_name` (str), the pipeline's name.

    Raises:
        ValueError: If `mode` is not "full" or "lite".
    """
    if mode not in ('full', 'lite'):
        raise ValueError(
            f"Invalid mode '{mode}'. Accepted values are 'full' or 'lite'.")

    config = {
        "path": {
            "config_dir":
            os.path.join(CBIG_REPDATA_DIR, 'stable_projects', 'preprocessing',
                         'Tan2026_TMSmotion', 'input_data', kwargs['pl_name'],
                         'preproc_rest', 'config'),
            "fmrinii_dir":
            os.path.join(CBIG_REPDATA_DIR, 'stable_projects', 'preprocessing',
                         'Tan2026_TMSmotion', 'input_data', kwargs['pl_name'],
                         'preproc_rest', 'fmrinii'),
            "rep_sel_subID_gt_dir":
            os.path.join(CBIG_REPDATA_DIR, 'stable_projects', 'preprocessing',
                         'Tan2026_TMSmotion', 'input_data', kwargs['pl_name'],
                         mode, 'select_subID_gt'),
            "rep_sel_subID_ngt_dir":
            os.path.join(CBIG_REPDATA_DIR, 'stable_projects', 'preprocessing',
                         'Tan2026_TMSmotion', 'input_data', kwargs['pl_name'],
                         mode, 'select_subID_ngt'),
            "rep_sel_subID_window_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'select_subID_window'),
            "rep_sel_subID_sim_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'select_subID_sim'),
            "rep_TMS_Kong2025_ngt_concat_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_ngt',
                         'concat'),
            "rep_TMS_Kong2025_ngt_outlier_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_ngt',
                         'outlier'),
            "rep_TMS_Kong2025_ngt_symlink_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_ngt',
                         'symlink'),
            "rep_TMS_Kong2025_ngt_TMS_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_ngt',
                         'TMS'),
            "rep_TMS_Kong2025_gt_gt_run_list_2_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_gt',
                         'gt_run_list_2'),
            "rep_TMS_Kong2025_gt_gt_split_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_gt',
                         'gt_split'),
            "rep_TMS_Kong2025_gt_TMS_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'TMS_Kong2025_gt',
                         'TMS'),
            "rep_compare_Kong2025_dep_dir":
            os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                         kwargs['pl_name'], 'output', 'compare_Kong2025_dep')
        },
        "pipeline": {
            "pl_name": kwargs['pl_name'],
            "mode": mode
        },
        "preproc": {
            "skip_N": kwargs['skip_N'],
            "FD_metric": kwargs['FD_metric'],
            "DV_metric": kwargs['DV_metric'],
            "FD_th": kwargs['FD_th'],
            "DV_th": kwargs['DV_th'],
            "discard_run": kwargs['discard_run'],
            "rm_seg": kwargs['rm_seg'],
            "bandstop": kwargs['bandstop'],
            "bandstop_TR_th": kwargs['bandstop_TR_th'],
            "bandstop_low_f": kwargs['bandstop_low_f'],
            "bandstop_high_f": kwargs['bandstop_high_f'],
            "lowpass_low_f": kwargs['lowpass_low_f'],
            "regressor": kwargs['regressor'],
            "bandpass_low_f": kwargs['bandpass_low_f'],
            "bandpass_high_f": kwargs['bandpass_high_f'],
            "FD_no_interp_th": kwargs['FD_no_interp_th'],
            "DV_no_interp_th": kwargs['DV_no_interp_th']
        },
        "select_subID": {
            "gt": {
                "min_gt_hr": kwargs['min_gt_hr'],
                "min_uniq_sesID": kwargs['min_uniq_sesID']
            },
            "window": {
                "window_time_sec":
                kwargs['window_time_sec'],
                "idx_frame_start_ini":
                kwargs['idx_frame_start_ini'],
                "percent_censored_bin_interval":
                kwargs['percent_censored_bin_interval'],
                "percent_censored_bin_lower":
                kwargs['percent_censored_bin_lower'],
                "percent_censored_bin_upper":
                kwargs['percent_censored_bin_upper'],
                "FD_metric_bin_interval":
                kwargs['FD_metric_bin_interval'],
                "FD_metric_bin_lower":
                kwargs['FD_metric_bin_lower'],
                "FD_metric_bin_upper":
                kwargs['FD_metric_bin_upper'],
                "DV_metric_bin_interval":
                kwargs['DV_metric_bin_interval'],
                "DV_metric_bin_lower":
                kwargs['DV_metric_bin_lower'],
                "DV_metric_bin_upper":
                kwargs['DV_metric_bin_upper'],
                "col_ngt_nframes_name":
                kwargs['col_ngt_nframes_name']
            },
            "sim": {
                "sim_iter":
                kwargs['sim_iter'],
                "num_window_per_run":
                kwargs['num_window_per_run'],
                "sample_bin_interval":
                kwargs['sample_bin_interval'],
                "sample_bin_min":
                kwargs['sample_bin_min'],
                "sample_bin_max":
                kwargs['sample_bin_max'],
                "num_sim_per_bin":
                kwargs['num_sim_per_bin'],
                "num_select_bin2select_subID":
                kwargs['num_select_bin2select_subID'],
                "num_sample_bin2select_subID":
                kwargs['num_sample_bin2select_subID'],
                "num_sim_run_per_sim_ses":
                kwargs['num_sim_run_per_sim_ses'],
                "sim_deficient_bin_iter":
                kwargs['sim_deficient_bin_iter'],
                "select_bin_min":
                kwargs['select_bin_min'],
                "select_bin_max":
                kwargs['select_bin_max'],
                "select_bin_interval":
                kwargs['select_bin_interval'],
                "rm_run_th":
                kwargs['rm_run_th'],
                "keep_select_bin_percent_th":
                kwargs['keep_select_bin_percent_th'],
                "pre_select_bin_interval":
                kwargs['pre_select_bin_interval'],
                "num_select_try":
                kwargs['num_select_try'],
                "lite_num_sim_ses":
                kwargs['lite_num_sim_ses'],
                "lite_combo":
                kwargs['lite_combo']
            }
        },
        "func2T1": {
            "func2T1_sm": kwargs['func2T1_sm']
        },
        "split": {
            "max_postcensored_min": kwargs['max_postcensored_min']
        },
        "TMS": {
            "FD_metric_value": kwargs['FD_metric_value'],
            "DV_metric_value": kwargs['DV_metric_value'],
            "oneBtwoA_flag": kwargs['oneBtwoA_flag'],
            "min5_flag": kwargs['min5_flag']
        },
        "MSHBM": {
            "MSHBM_w": kwargs['MSHBM_w'],
            "MSHBM_MRF": kwargs['MSHBM_MRF']
        },
        "compare": {
            "compare_num_sim_run_per_sim_ses":
            kwargs['compare_num_sim_run_per_sim_ses']
        },
        "stem": {
            "FD_stem": kwargs['FD_stem'],
            "DV_stem": kwargs['DV_stem']
        }
    }

    # Write config to JSON file
    output_dir = os.path.dirname(
        f"{TAN2026_TMSMOTION_DIR}/replication/{kwargs['pl_name']}/config")
    os.makedirs(output_dir, exist_ok=True)
    output_filename = f"{TAN2026_TMSMOTION_DIR}/replication/{kwargs['pl_name']}/config/pipeline.json"
    util.write_config_dict2json(config, output_filename)


def read_pl_ini(ini_file):
    """Read a pipeline `.ini` configuration file into a parameter dictionary.

    Parses all pipeline parameters from the `[pipeline]`, `[preproc]`,
    `[select_subID]`, `[func2T1]`, `[split]`, `[TMS]`, `[MSHBM]`,
    `[compare]`, `[stem]`, and optional `[lite_combo]` sections, casting
    each value to its expected type (int/float/str/bool/list). Option
    name case is preserved (not lowercased) so that subID-like keys
    under `[lite_combo]` (e.g. `sub-MSC07`) round-trip correctly.

    Args:
        ini_file (str): Path to the `.ini` file containing pipeline
            configuration.

    Returns:
        dict: A dictionary with all parameters extracted from the `.ini`
            file, keyed by parameter name (e.g. `pl_name`, `FD_metric`,
            `rm_run_th`).
    """
    config = configparser.ConfigParser()
    # Preserve option name case (needed for `[lite_combo]` subIDs like
    # `sub-MSC07`, which default ConfigParser would otherwise lowercase)
    config.optionxform = str
    config.read(ini_file)

    # Use dictionary to dynamically retrieve all parameters
    pl_ini_dict = {
        "pl_name":
        str(config['pipeline']['pl_name']),
        "skip_N":
        int(config['preproc']['skip_N']),
        "FD_metric":
        str(config['preproc']['FD_metric']),
        "DV_metric":
        str(config['preproc']['DV_metric']),
        "FD_th":
        float(config['preproc']['FD_th']),
        "DV_th":
        int(config['preproc']['DV_th']),
        "discard_run":
        int(config['preproc']['discard_run']),
        "rm_seg":
        int(config['preproc']['rm_seg']),
        "bandstop":
        config.getboolean('preproc', 'bandstop'),
        "bandstop_TR_th":
        float(config['preproc']['bandstop_TR_th']),
        "bandstop_low_f":
        float(config['preproc']['bandstop_low_f']),
        "bandstop_high_f":
        float(config['preproc']['bandstop_high_f']),
        "lowpass_low_f":
        float(config['preproc']['lowpass_low_f']),
        "regressor":
        str(config['preproc']['regressor']),
        "bandpass_low_f":
        float(config['preproc']['bandpass_low_f']),
        "bandpass_high_f":
        float(config['preproc']['bandpass_high_f']),
        "FD_no_interp_th":
        int(config['preproc']['FD_no_interp_th']),
        "DV_no_interp_th":
        int(config['preproc']['DV_no_interp_th']),
        "sim_iter":
        int(float(config['select_subID']['sim_iter'])),
        "min_gt_hr":
        int(config['select_subID']['min_gt_hr']),
        "min_uniq_sesID":
        int(config['select_subID']['min_uniq_sesID']),
        "num_window_per_run":
        int(config['select_subID']['num_window_per_run']),
        "sample_bin_interval":
        int(config['select_subID']['sample_bin_interval']),
        "sample_bin_min":
        int(config['select_subID']['sample_bin_min']),
        "sample_bin_max":
        int(config['select_subID']['sample_bin_max']),
        "num_sim_per_bin":
        int(config['select_subID']['num_sim_per_bin']),
        "lite_num_sim_ses":
        int(config['select_subID']['lite_num_sim_ses']),
        "lite_combo":
        dict(config['lite_combo']) if config.has_section('lite_combo') else {},
        "num_select_bin2select_subID":
        int(config['select_subID']['num_select_bin2select_subID']),
        "num_sample_bin2select_subID":
        int(config['select_subID']['num_sample_bin2select_subID']),
        "num_sim_run_per_sim_ses":
        int(config['select_subID']['num_sim_run_per_sim_ses']),
        "sim_deficient_bin_iter":
        int(float(config['select_subID']['sim_deficient_bin_iter'])),
        "select_bin_min":
        int(config['select_subID']['select_bin_min']),
        "select_bin_max":
        int(config['select_subID']['select_bin_max']),
        "select_bin_interval":
        int(config['select_subID']['select_bin_interval']),
        "window_time_sec":
        int(config['select_subID']['window_time_sec']),
        "idx_frame_start_ini":
        int(config['select_subID']['idx_frame_start_ini']),
        "percent_censored_bin_interval":
        float(config['select_subID']['percent_censored_bin_interval']),
        "percent_censored_bin_lower":
        float(config['select_subID']['percent_censored_bin_lower']),
        "percent_censored_bin_upper":
        float(config['select_subID']['percent_censored_bin_upper']),
        "FD_metric_bin_interval":
        float(config['select_subID']['FD_metric_bin_interval']),
        "FD_metric_bin_lower":
        float(config['select_subID']['FD_metric_bin_lower']),
        "FD_metric_bin_upper":
        float(config['select_subID']['FD_metric_bin_upper']),
        "DV_metric_bin_interval":
        float(config['select_subID']['DV_metric_bin_interval']),
        "DV_metric_bin_lower":
        float(config['select_subID']['DV_metric_bin_lower']),
        "DV_metric_bin_upper":
        float(config['select_subID']['DV_metric_bin_upper']),
        "col_ngt_nframes_name":
        str(config['select_subID']['col_ngt_nframes_name']),
        "rm_run_th":
        list(map(int, config['select_subID']['rm_run_th'].split(','))),
        "keep_select_bin_percent_th":
        int(config['select_subID']['keep_select_bin_percent_th']),
        "pre_select_bin_interval":
        int(config['select_subID']['pre_select_bin_interval']),
        "num_select_try":
        int(config['select_subID']['num_select_try']),
        "func2T1_sm":
        int(config['func2T1']['func2T1_sm']),
        "max_postcensored_min":
        int(config['split']['max_postcensored_min']),
        "FD_metric_value":
        list(map(float, config['TMS']['FD_metric_value'].split(','))),
        "DV_metric_value":
        list(map(int, config['TMS']['DV_metric_value'].split(','))),
        "oneBtwoA_flag":
        str(config['TMS']['oneBtwoA_flag']),
        "min5_flag":
        str(config['TMS']['min5_flag']),
        "MSHBM_w":
        int(config['MSHBM']['MSHBM_w']),
        "MSHBM_MRF":
        int(config['MSHBM']['MSHBM_MRF']),
        "compare_num_sim_run_per_sim_ses":
        int(config['compare']['compare_num_sim_run_per_sim_ses']),
        "FD_stem":
        str(config['stem']['FD_stem']),
        "DV_stem":
        str(config['stem']['DV_stem'])
    }

    # Return all parameters as unpacked dictionary
    return pl_ini_dict


def main():
    """Generate a pipeline configuration JSON file from a `.ini` input file.

    Reads the main config, then the pipeline `.ini` file and mode given
    as the first and second command-line arguments (`sys.argv[1]`,
    `sys.argv[2]`), and uses them to generate and write the pipeline's
    `pipeline.json` configuration via `generate_pl_config`.
    """
    # Read main config
    main_config = os.path.join(TAN2026_TMSMOTION_DIR, 'replication', 'config',
                               'main_config.json')
    util.read_json_file(main_config)

    # Read user-defined input
    ini_file = sys.argv[1]
    mode = sys.argv[2]

    # Read .ini as input to generate_pl_config
    pl_ini_dict = read_pl_ini(ini_file)

    # Generate the pipeline configuration
    generate_pl_config(TAN2026_TMSMOTION_DIR, mode, **pl_ini_dict)


if __name__ == "__main__":
    main()

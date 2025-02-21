"""
This script modifies several configuration and setup scripts within the CBIG codebase.
It defines functions to replace specific code blocks in various files with new code blocks or remove them entirely.

Written by Tian Fang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
import os
import re


def replace_code(file_path, code_start, code_end, new_code):
    """
    Replace a code block in a file with a new code block.

    Parameters
    ----------
    file_path : str
        Path to the file to modify.
    code_start : str
        Start of the code block to replace.
    code_end : str
        End of the code block to replace.
    new_code : str
        New code block to insert.
    """
    with open(file_path, 'r') as file:
        content = file.read()

    # Use regular expressions to replace the code block
    escaped_start = re.escape(code_start)
    escaped_end = re.escape(code_end)

    # check whether escaped code_start is found
    if not re.search(escaped_start, content):
        print("Code start not found")
        return

    # check whether escaped code_end is found
    if not re.search(escaped_end, content):
        print("Code end not found")
        return

    pattern = re.compile(escaped_start + r'.*?' + escaped_end, re.DOTALL)
    modified_content = re.sub(pattern, new_code, content)

    # Check if changes were made
    if content == modified_content:
        print("No changes were made.")
        return

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.write(modified_content)

    print("Code block replaced successfully.")


def modify_generic_setup(cbig_code_dir):
    file_path = os.path.join(cbig_code_dir, 'setup', 'CBIG_generic_setup.sh')

    code_start = '############################\n# set up Workbench directory'
    code_end = 'to check if the results can be replicated."\n  fi\nfi'
    new_code = ''
    replace_code(file_path, code_start, code_end, new_code)

    code_start = '####################################\n# create symlink to git hook scripts'
    code_end = 'ln -s "$CBIG_CODE_DIR/hooks/pre-commit" "$CBIG_CODE_DIR/.git/hooks/pre-commit"\nfi'
    new_code = ''
    replace_code(file_path, code_start, code_end, new_code)

    code_start = "# check whether user's default matlab"
    code_end = 'to check if the results can be replicated."\n  fi'
    new_code = ''
    replace_code(file_path, code_start, code_end, new_code)


def modify_pbsubmit(cbig_code_dir):
    file_path = os.path.join(cbig_code_dir, 'setup', 'CBIG_pbsubmit')
    code_start = '# Check current node'
    code_end = 'sleep ${sleeptime}'
    new_code = 'cd ${work_dir}\neval "${cmd_script}"'
    replace_code(file_path, code_start, code_end, new_code)


def modify_main_script(cbig_code_dir):
    file_path = os.path.join(cbig_code_dir, 'stable_projects', 'preprocessing',
                             'CBIG_fMRI_Preproc2016',
                             'CBIG_preproc_fMRI_preprocess.csh')
    code_start = 'set matlab_runtime_util'
    code_end = '# MATLAB Runtime utilities folder'
    new_code = (
        'set matlab_runtime_util = '
        # ! Replace the following line if the compiled utilities is not under the container directory
        '"$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/utilities" '
        '# MATLAB Runtime utilities folder')
    replace_code(file_path, code_start, code_end, new_code)


def modify_single_unit_test(cbig_code_dir):
    file_path = os.path.join(cbig_code_dir, 'stable_projects', 'preprocessing',
                             'CBIG_fMRI_Preproc2016', 'unit_tests',
                             'single_subject',
                             'CBIG_preproc_single_subject_unit_test.m')
    code_start = '%% we need to periodically check whether the job has finished or not'
    code_end = 'pause(60); % sleep for 1min and check again\n            end'
    new_code = ''
    replace_code(file_path, code_start, code_end, new_code)


def modify_config(cbig_code_dir):
    file_path = os.path.join('app', 'setup', 'CBIG_config_in_container.sh')
    code_start = 'export CBIG_CODE_DIR'
    code_end = 'storage/CBIG_private'
    new_code = 'export CBIG_CODE_DIR=/app/CBIG_private'
    replace_code(file_path, code_start, code_end, new_code)

    # Then append a block of code to the end of the file
    with open(file_path, 'a') as file:
        appended_code = """

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/app/miniconda/bin/conda' 'shell.bash' 'hook' 2>/dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/app/miniconda/etc/profile.d/conda.sh" ]; then
        . "/app/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/app/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

export PATH="/app/miniconda/bin:$PATH"

conda activate CBIG_py3

export PREPROC2016=${CBIG_CODE_DIR}/stable_projects/preprocessing
export PREPROC2016=${PREPROC2016}/CBIG_fMRI_Preproc2016/CBIG_preproc_fMRI_preprocess.csh

############################################
# NEW: Set up MATLAB Compiler Runtime
############################################
export APPS_DIR=/apps
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/runtime/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/bin/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/sys/os/glnxa64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${APPS_DIR}/matlab/MCR_R2018b/v95/sys/opengl/lib/glnxa64
export MATLAB_RUNTIME_DIR=${APPS_DIR}/matlab/MCR_R2018b/v95

############################################
# NEW: Add ${CBIG_CODE_DIR}/utilities/scripts/ to PATH
############################################
export PATH=${CBIG_CODE_DIR}/utilities/scripts/:$PATH

"""
        file.write(appended_code)


if __name__ == "__main__":
    # Get the environment variable
    cbig_code_dir = os.getenv('CBIG_CODE_DIR', '/app/CBIG_private')
    print(f'CBIG_CODE_DIR: {cbig_code_dir}')

    # Modify existing scripts in-place
    modify_generic_setup(cbig_code_dir)
    modify_pbsubmit(cbig_code_dir)
    modify_main_script(cbig_code_dir)
    modify_single_unit_test(cbig_code_dir)
    modify_config(cbig_code_dir)

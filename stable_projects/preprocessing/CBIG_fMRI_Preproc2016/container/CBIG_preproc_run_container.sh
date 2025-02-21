#!/bin/bash

# This script is a wrapper for running the CBIG fMRI preprocessing pipeline
# with containerization technology, namely docker and singularity.
# It performs the following tasks:
# - Parses command-line arguments to determine container type (Docker or Singularity),
#   paths for MATLAB Compiler Runtime (MCR), input, and output directories.
# - Creates the output directory if it does not exist.
# - Pulls the Docker image from Docker Hub and converts it to a Singularity image if necessary.
# - Constructs and executes the appropriate Docker or Singularity command
#   to run the preprocessing pipeline with the provided arguments.
# Written by Tian Fang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Function to print usage
usage() {
    echo -e "Usage: $0 [WRAPPER_ARG]...[MAIN_SCRIPT_ARG]... \n\n\
This CBIG fMRI Preprocessing wrapper script does the following things: \n\
- creates the output directory if it does not exist \n\
- pulls the docker image from Docker Hub (when using Docker) \n\
- converts the Docker image to Singularity image (when using Singularity) \n\
- generates Docker/Singularity commands for you \n\
- executes the Docker/Singularity commands and passes [MAIN_SCRIPT_ARG] to the main preprocessing script\n"
    echo "  -c, --container        Containerization technology: 'docker' or 'singularity'"
    echo "  -m, --mcr              Path to MATLAB Compiler Runtime (MCR)"
    echo "  -i, --input            Path to input directory"
    echo "  -o, --out              Path to output directory"
    echo "  -sdir, --singularity-dir  Optional: Directory to store Singularity image when using Singularity,\
default is the current directory"
    echo "  -h, --help             Optional: Display this help and exit"
    echo "  -mh, --main-help       Optional: Display help for the main preprocessing script"
    echo "  [MAIN_SCRIPT_ARG] additional arguments other than the above \
will be passed directly to the main preprocessing script in the container."
    echo -e "\nExample:\n\
${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/CBIG_preproc_run_container.sh \\
--container singularity \\
--mcr /apps/matlab/MATLAB_Compiler_Runtime/v95 \\
--input \$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject \\
--out \${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/tmp/singularity_out \\
-sdir \${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/tmp/ \\
-s sub005 \\
-anat_s sub005_FS \\
-anat_d /input/data \\
-fmrinii /input/docker/sub005_docker.fmrinii \\
-config /input/docker/prepro_docker.config \\
-nocleanup"
    exit 1
}

pull_docker() {
    echo "Pulling the Docker image thomasyeolab/cbig_fmri_preproc2016:latest..."
    docker_pull_cmd="docker pull thomasyeolab/cbig_fmri_preproc2016:latest"
    echo "Docker pull command: $docker_pull_cmd"
    echo -e "\n\n"
    eval $docker_pull_cmd
}

run_docker() {
    path_to_MCR=$1
    path_to_input=$2
    path_to_output=$3
    shift 3
    main_script_arg=$@

    echo "Running the Docker container..."
    docker_cmd="docker run -it --rm \
        -u $(id -u):$(id -g) -e HOME=$HOME \
        -v \"$path_to_MCR\":/apps/matlab/MCR_R2018b/v95:ro \
        -v \"$path_to_input\":/input:ro \
        -v \"$path_to_output\":/out \
        thomasyeolab/cbig_fmri_preproc2016:latest \
        -output_d /out $main_script_arg"
    echo "Docker command: $docker_cmd"
    echo -e "\n\n"
    eval $docker_cmd
}

invalid_container() {
    echo "Invalid containerization technology selected. Please choose 'docker' or 'singularity'."
    usage
}

build_singularity() {
    singularity_image=$1

    # Convert Docker image to Singularity image if it does not exist
    if [[ ! -d $singularity_image ]]; then
        # Create Singularity image directory if it does not exist
        if [[ ! -d $(dirname $singularity_image) ]]; then
            echo "Creating Singularity image directory: $(dirname $singularity_image)"
            mkdir -p $(dirname $singularity_image)
        fi
        echo "Converting Docker image to Singularity image (takes about half an hour)..."
        singularity_cmd="singularity build --sandbox \"$singularity_image\" \
docker://thomasyeolab/cbig_fmri_preproc2016:latest"
        echo -e "Singularity build command: $singularity_cmd"
        echo -e "\n\n"
        eval $singularity_cmd
    fi
}

run_singularity() {
    path_to_MCR=$1
    path_to_input=$2
    path_to_output=$3
    singularity_image=$4
    shift 4
    main_script_arg=$@

    echo "Running the Singularity container..."
    singularity_run_cmd="singularity run --cleanenv --no-home --writable \
        --bind \"$path_to_MCR\":/apps/matlab/MCR_R2018b/v95:ro \
        --bind \"$path_to_input\":/input:ro \
        --bind \"$path_to_output\":/out \
        \"$singularity_image\" \
        -output_d /out $main_script_arg"
    echo "Singularity run command: $singularity_run_cmd"
    echo -e "\n\n"
    eval $singularity_run_cmd
}
# Parse command line arguments
singularity_dir=""
main_script_arg=""
display_main_help=false

while [[ "$#" -gt 0 ]]; do
    case $1 in
    -c | --container)
        container_type="$2"
        shift
        ;;
    -m | --mcr)
        path_to_MCR="$2"
        shift
        ;;
    -i | --input)
        path_to_input="$2"
        shift
        ;;
    -o | --out)
        path_to_output="$2"
        shift
        ;;
    -sdir | --singularity-dir)
        singularity_dir="$2"
        shift
        ;;
    -h | --help)
        usage
        ;;
    -mh | --main-help)
        display_main_help=true
        ;;
    *)
        main_script_arg+="$1 $2 "
        shift
        ;;
    esac
    shift
done

# Display help for the main preprocessing script if requested
if [[ $display_main_help == true ]]; then
    # if container_type is empty, check whether docker or singularity is installed
    if [[ -z $container_type ]]; then
        if command -v docker &>/dev/null; then
            container_type="docker"
        elif command -v singularity &>/dev/null; then
            container_type="singularity"
        else
            echo "Docker or Singularity is required to run the fMRI preprocessing pipeline."
            exit 1
        fi
    fi
    if [[ $container_type == "docker" ]]; then
        pull_docker
        docker run -it --rm thomasyeolab/cbig_fmri_preproc2016:latest
    elif [[ $container_type == "singularity" ]]; then
        # Define Singularity image path
        singularity_image="${singularity_dir:-$PWD}/cbig_fmri_preproc2016"

        build_singularity $singularity_image
        singularity run --cleanenv $singularity_image
    else
        invalid_container
    fi
    exit 0
fi

# Check if required arguments are provided
if [[ -z $container_type || -z $path_to_MCR || -z $path_to_input || -z $path_to_output ]]; then
    usage
fi

# Create output directory if it does not exist
if [[ ! -d $path_to_output ]]; then
    echo "Creating output directory: $path_to_output"
    mkdir -p "$path_to_output"
fi

# Run the appropriate container command
if [[ $container_type == "docker" ]]; then
    docker_pull_cmd
    run_docker $path_to_MCR $path_to_input $path_to_output $main_script_arg
elif [[ $container_type == "singularity" ]]; then
    # Define Singularity image path
    singularity_image="${singularity_dir:-$PWD}/cbig_fmri_preproc2016"

    build_singularity $singularity_image
    run_singularity $path_to_MCR $path_to_input $path_to_output $singularity_image $main_script_arg
else
    invalid_container
fi

#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -d <docs> -t <setting> -k <no_topics> -r <run_file> -m <polarLDA_dir> -o <out_dir> [-q <queue>]
    - docs         Text file with each line summarizing the document of first 
                    modality; e.g., brain atrophy doc
    - setting       Settings file for polarLDA model; e.g., 
                    CBIG_ASDf_polarLDA_infSettings.txt
    - no_topics     Number of topics / factors; e.g., 3
    - run_file      a file containing a run number in each row; e.g., 1, 2, 3, ..., 100
    - polarLDA_dir  Directory of binary executable file for polarLDA 
    - out_dir       Output directory; e.g., ~/outputs/polarLDA/
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":d:t:k:r:m:o:q:" opt; do
    case "${opt}" in
        d) docs=${OPTARG};;
        t) setting=${OPTARG};;
        k) no_topics=${OPTARG};;
        r) run_file=${OPTARG};;
        m) polarLDA_dir=${OPTARG};;
        o) out_dir=${OPTARG};;
        q) queue=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${docs}" ] || [ -z "${setting}" ] || \
     [ -z "${no_topics}" ] || [ -z "${run_file}" ] || \
     [ -z "${polarLDA_dir}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

echo '---polarLDA estimation.'

mkdir -p ${out_dir}/k${no_topics}

progress_file=${out_dir}/k${no_topics}/progress.txt
> ${progress_file}
for r in `cat ${run_file}`; do
    run_dir=${out_dir}/k${no_topics}/r${r}
    rm -rf ${run_dir}
    mkdir -p ${run_dir}
    log_file=${run_dir}/polarLDA.log
    > ${log_file}

    # initialize alpha to be 1/no_topics
    alpha=$(echo 1/${no_topics} | bc -l)

    if [ -z "${queue}" ]; then
        # converting relative paths to absolute for qsub
        setting=$(readlink -f ${setting})
        docs=$(readlink -f ${docs})
        run_dir=$(readlink -f ${run_dir})

        date >> ${log_file}
        echo "Docs: ${docs} " >> ${log_file}
        echo "Number of topics: ${no_topics}" >> ${log_file}
        echo "Settings:" >> ${log_file}
        cat ${setting} >> ${log_file}

        ${polarLDA_dir}/polarLDA est ${alpha} ${no_topics} ${setting} ${docs} random ${run_dir}/ ${r} >> ${log_file}
        echo "${r}" >> ${progress_file}
    else
    
        # converting relative paths to absolute for qsub
        setting=$(readlink -f ${setting})
        docs=$(readlink -f ${docs})
        run_dir=$(readlink -f ${run_dir})

        date >> ${log_file}
        echo "Docs: ${docs} " >> ${log_file}
        echo "Number of topics: ${no_topics}" >> ${log_file}
        echo "Settings:" >> ${log_file}
        cat ${setting} >> ${log_file}
        cmd="${polarLDA_dir}/polarLDA est ${alpha} ${no_topics} ${setting} ${docs} \
        random ${run_dir}/ ${r} >> ${log_file}"
        cmd="${cmd};echo "${r}" >> ${progress_file}"
        # write command to a script file
        script_file=${run_dir}/job_submitted.sh
        echo "#!/bin/bash" > ${script_file}
        echo $cmd >> ${script_file}
        chmod 750 ${script_file}
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${script_file}" -walltime 10:00:00 -mem 16G -name "polarLDA_est" \
        -joberr ${run_dir}/polarLDA.err -jobout ${run_dir}/polarLDA.out
    fi
done


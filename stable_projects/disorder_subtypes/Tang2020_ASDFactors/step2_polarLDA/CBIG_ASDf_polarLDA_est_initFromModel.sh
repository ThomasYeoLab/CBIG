#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -d <docs> -t <setting> -k <no_topics> -i <init> -r <run_file> \
-m <model_name> -p <polarLDA_dir> -o <out_dir> [-q <queue>]
    - docs         Text file with each line summarizing the document;
    - setting       Setting file for polarLDA model; e.g., 
                    ./CBIG_ASDf_polarLDA_infSettings.txt
    - no_topics     Number of topics / factors; e.g., 3
    - init          'random' or 'model'. 
                    'random' means random initialization. run_file argument is required.
                    'model' means initialization with given model parameters. run_file argument is not required.
    - run_file      (Optional) a file containing a run number in each row; e.g., 1, 2, 3, ..., 100
    - model_name    (Optional) absolute path to a model; e.g., ~/output/polarLDA/k3/final
    - polarLDA_dir  Directory of binary executable file for polarLDA 
    - out_dir       Output directory; e.g., ~/outputs/polarLDA/
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":d:t:k:i:r:m:p:o:q:" opt; do
    case "${opt}" in
        d) docs=${OPTARG};;
        t) setting=${OPTARG};;
        k) no_topics=${OPTARG};;
        i) init=${OPTARG};;
        r) run_file=${OPTARG};;
        m) model_name=${OPTARG};;
        p) polarLDA_dir=${OPTARG};;
        o) out_dir=${OPTARG};;
        q) queue=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${docs}" ] || [ -z "${setting}" ] || \
     [ -z "${no_topics}" ] || [ -z "${init}" ] || \
     [ -z "${polarLDA_dir}" ] || [ -z "${out_dir}" ]; then
    echo Missing Parameters!
    usage
fi
if [ "${init}" == "random" ] && [ -z "${run_file}" ]; then
    echo Missing run_file argument.
    usage
fi
if [ "${init}" == "model" ] && [ -z "${model_name}" ]; then
    echo Missing model_name argument.
    usage
fi

###########################################
# Main
###########################################

echo '---polarLDA estimation.'

mkdir -p ${out_dir}/k${no_topics}

progress_file=${out_dir}/k${no_topics}/progress.txt
> ${progress_file}

if [ "${init}" == "random" ]; then
    for r in `cat ${run_file}`; do
        run_dir=${out_dir}/k${no_topics}/r${r}
        rm -r ${run_dir}
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
            qsub -V -q ${queue} << EOJ
#!/bin/bash
#PBS -N 'polarLDA_est'
#PBS -P 11000481
#PBS -l walltime=60:00:0
#PBS -l mem=16gb
#PBS -e ${run_dir}/polarLDA.err
#PBS -o ${run_dir}/polarLDA.out
    
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
EOJ
        fi
    done
elif [ "${init}" == "model" ]; then
    run_dir=${out_dir}/k${no_topics}/init_with_model
    rm -r ${run_dir}
    mkdir -p ${run_dir}
    log_file=${run_dir}/polarLDA.log
    > ${log_file}

    # initialize alpha to be 1/no_topics
    # this is a dummy alpha, the alpha will be initialized with model alpha
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

        ${polarLDA_dir}/polarLDA est ${alpha} ${no_topics} ${setting} ${docs} ${model_name} ${run_dir}/ >> ${log_file}
        echo "Done" >> ${progress_file}
    else
        qsub -V -q ${queue} << EOJ
#!/bin/bash
#PBS -N 'polarLDA_est'
#PBS -l walltime=60:00:0
#PBS -l mem=16gb
#PBS -e ${run_dir}/polarLDA.err
#PBS -o ${run_dir}/polarLDA.out
    
    # converting relative paths to absolute for qsub
    setting=$(readlink -f ${setting})
    docs=$(readlink -f ${docs})
    run_dir=$(readlink -f ${run_dir})

    date >> ${log_file}
    echo "Docs: ${docs} " >> ${log_file}
    echo "Number of topics: ${no_topics}" >> ${log_file}
    echo "Settings:" >> ${log_file}
    cat ${setting} >> ${log_file}

    ${polarLDA_dir}/polarLDA est ${alpha} ${no_topics} ${setting} ${docs} ${model_name} ${run_dir}/ >> ${log_file}
    echo "Done" >> ${progress_file}
EOJ
    fi
fi


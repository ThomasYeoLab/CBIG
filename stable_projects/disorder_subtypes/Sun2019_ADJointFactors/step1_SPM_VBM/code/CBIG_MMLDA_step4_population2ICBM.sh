#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <inDir> -s <scriptDir> -p <spmDir> [-q <queue>]
    Population the dartel template to ICBM space.
    - inDir         Input directory of T1 images
    - scriptDir     Script directory of current script
    - spmDir        SPM installation directory
    - queue         (Optional) if you have a cluster, use it to specify the 
                    queue to which you want to qsub these jobs; if not provided,
                    jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:s:p:q:" opt; do
    case "${opt}" in
            i) inDir=${OPTARG};;
            s) scriptDir=${OPTARG};;
            p) spmDir=${OPTARG};;
            q) queue=${OPTARG};;
            *) usage;;
        esac
done
shift $((OPTIND-1))
if [ -z "${inDir}" ] || [ -z "${scriptDir}" ] || [ -z "${spmDir}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################
echo 'Step 4: Population to ICBM.'

mkdir -p ${inDir}
logDir=${inDir}/logs/step4_population2ICBM
mkdir -p ${logDir}
progressFile=${logDir}/progress.txt
> ${progressFile}

id=population2ICBM
logFile=${logDir}/${id}.log
imgPath=${inDir}/${id}.nii
if [ -z "${queue}" ]; then
    matlab -nodisplay -nosplash -r \
    "spm_dir='${spmDir}';script_dir='${scriptDir}';\
    output_dir='${inDir}';CBIG_MMLDA_population_to_ICBM;exit;" > ${logFile}
    echo "Done" >> ${progressFile}
else
    qsub -V -q ${queue} << EOJ

#!/bin/sh
#PBS -N 'step4_population2ICBM'
#PBS -l walltime=1:00:0
#PBS -l mem=8gb   
#PBS -e ${logDir}/${id}.err
#PBS -o ${logDir}/${id}.out
    
    cd ${scriptDir}
    
    matlab -nodisplay -nosplash -r \
    "spm_dir='${spmDir}';script_dir='${scriptDir}';\
    output_dir='${inDir}';CBIG_MMLDA_population_to_ICBM;exit;" > ${logFile}
    echo "Done" >> ${progressFile}
EOJ
fi

./CBIG_MMLDA_waitUntilFinished.sh ${progressFile} 

echo 'Step 4: Population to ICBM. -- Finished.'

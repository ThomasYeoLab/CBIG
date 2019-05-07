#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

progress_file=${1}
total_no_jobs=${2}

# Do not return until all jobs are finished
no_jobs_done=0
no_jobs_done_prev=-1
while [ ${no_jobs_done} -lt ${total_no_jobs} ]; do
    if [ ! ${no_jobs_done} -eq ${no_jobs_done_prev} ]; then
        echo "${no_jobs_done}/${total_no_jobs} Finished"
    fi
    sleep 5
    no_jobs_done_prev=${no_jobs_done}
    no_jobs_done=`grep -c "^" ${progress_file}`
done
#!/bin/bash

# Wait jobs with <prefix> to finish
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

shopt -s expand_aliases
sleep_time=1m
me=$(whoami)
alias myqstat="qstat -u $me | grep "${1}" | wc -l"
status=$(myqstat)
while [ $status -ne 0 ]; do
    echo "Waitting for existing jobs with prefix <"${1}"> to finish..."
    sleep $sleep_time
    status=$(myqstat)
done

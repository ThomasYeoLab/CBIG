#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

shopt -s expand_aliases
me=$(whoami)
#alias myqstat="qstat -u $me | grep "${1}" | wc -l"
alias myqstat="qstat -f | grep -C 1 ${me} | grep ${1} | wc -l"
status=$(myqstat)
if [ $status -ne 0 ]; then
	echo "Waitting for existing "${1}" jobs to finish..."
fi
while [ $status -ne 0 ]; do
    sleep 10s
    status=$(myqstat)
done

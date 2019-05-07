#!/bin/sh

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

progressFile=${1}
fullList=${2}

# Do not return until all jobs are finished
if [ -z ${fullList} ]; then 
    noJobs=1
else
    noJobs=`grep -c "^" ${fullList}`
fi
noJobsDone=0
noJobsDonePrev=-1
while [ ${noJobsDone} -lt ${noJobs} ]; do
	if [ ! ${noJobsDone} -eq ${noJobsDonePrev} ]; then
		echo ${noJobsDone}/${noJobs} Finished
	fi
	sleep 5
	noJobsDonePrev=${noJobsDone}
	noJobsDone=`grep -c "^" ${progressFile}`
done

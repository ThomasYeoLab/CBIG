#!/bin/sh

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

queue=circ-spool

# BET on baseline scans in 11 batches
# Each batch uses a different set of BET parameters/procedure
outDir=~/storage/forPNASRelease/outputs/VBM_bl/brains/
for i in {1..11}; do
	imgList=./imgList_BETParam_bl/imgList${i}.txt
	BETparam=./imgList_BETParam_bl/CBIG_BETParam${i}.sh
	../lib/CBIG_runBET.sh -i ${imgList} -p ${BETparam} -o ${outDir} -q ${queue}
done


# BET on m24 scans in 10 batches
# Each batch uses a different set of BET parameters/procedure
outDir=~/storage/forPNASRelease/outputs/VBM_m24/brains/
for i in {1..10}; do
	imgList=./imgList_BETParam_m24/imgList${i}.txt;
	BETparam=./imgList_BETParam_m24/CBIG_BETParam${i}.sh
	../lib/CBIG_runBET.sh -i ${imgList} -p ${BETparam} -o ${outDir} -q ${queue}
done

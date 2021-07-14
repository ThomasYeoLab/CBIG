#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

file=$1
num=$2
nfinished=`cat $file | wc -l`
while [ $nfinished -lt ${num} ]
do
nfinished=`cat $file | wc -l`
njob=`qstat | grep $(whoami) | wc -l`
if [ "$njob" -eq 0 ];then
	break
fi
sleep 1m
done

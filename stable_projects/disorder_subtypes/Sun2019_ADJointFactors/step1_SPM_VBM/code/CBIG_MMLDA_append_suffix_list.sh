#!/bin/bash

# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

inList=${1}
suffix=${2}
outList=${3}

# append suffix to each line of a list
sed "s/$/${suffix}/" ${inList} > ${outList}


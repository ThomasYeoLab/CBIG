#!/bin/bash
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=`readlink -f $0 | xargs dirname`/../
python -m unittest discover --verbose --start-directory $ROOTDIR/cbig/Nguyen2020/

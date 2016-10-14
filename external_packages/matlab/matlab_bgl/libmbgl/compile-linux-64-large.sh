#!/bin/bash -e

BOOST_DIR=${HOME}/dev/lib/boost_1_36_0/
YASMIC_DIR=.

source ccfiles.sh
OFILES=`echo ${CCFILES} | sed -e 's/\.cc/\.o/g'`

CFLAGS="-O2 -DMATLAB_BGL_LARGE_ARRAYS -fPIC -c -I${BOOST_DIR} -I${YASMIC_DIR}"
#CFLAGS="-g -W -DMATLAB_BGL_LARGE_ARRAYS -fPIC -c -I${BOOST_DIR} -I${YASMIC_DIR}"

function echocmd {
    echo $@
    $@
}

for file in ${CCFILES}; do
    echocmd g++-3.4 $CFLAGS $file
done

echocmd ar rc libmbgl-linux-64-large.a ${OFILES} 
    
echocmd rm ${OFILES}    

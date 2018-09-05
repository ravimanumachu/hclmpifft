#!/bin/bash

###################################################
FFTWRESDIR=fftw2dresults_t72
mkdir -p ${FFTWRESDIR}

for (( COUNT=16; COUNT<=1000; COUNT++ ));
do
    N=`expr $COUNT \* 64`

    ./fftw $N $N 72 0 >> ${FFTWRESDIR}/fftw3_t72.res
    sleep 2

done

###################################################

exit 0

###################################################

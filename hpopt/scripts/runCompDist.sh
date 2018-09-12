#!/bin/bash

#################################################################

BLASBASEDIR=Results/
DGEMMETIMEDIR=${BLASBASEDIR}/parts
mkdir -p ${DGEMMETIMEDIR}

#################################################################

for (( COUNT=1; COUNT<=2100; COUNT+=1 ));
do
        N=`expr $COUNT \* 64`

				#echo ./compDist $N 1024 0 /home/hamidreza/publications/2017/heterogeneous-energy-performance/Files/funcntions/dgemm/compare/

        ./compDist $N 4736 1 1 /home/hamidreza/publications/2017/heterogeneous-performance/Files/funcntions/dgemm_sep/compare/ > ${DGEMMETIMEDIR}/${N} 2>&1
        
        TIMETH=`cat ${DGEMMETIMEDIR}/${N} | grep "Time Threshold (Original)=" | awk '{print $4}'`
        HETTIME=`cat ${DGEMMETIMEDIR}/${N} | grep "Heterogeneous_Original-> ExecTime=" | awk '{print $3}'`
        USEDPROC_HET=`cat ${DGEMMETIMEDIR}/${N} | grep "Heterogeneous_Original->" | awk '{print $9}'`
        
        SMOOTHTIME=`cat ${DGEMMETIMEDIR}/${N} | grep "Heterogeneous_Smooth-> ExecTime=" | awk '{print $3}'`
        USEDPROC_SMOOTH=`cat ${DGEMMETIMEDIR}/${N} | grep "Heterogeneous_Smooth->" | awk '{print $9}'`
        
        CONSTTIME=`cat ${DGEMMETIMEDIR}/${N} | grep "ConstantPM-> ExecTime=" | awk '{print $3}'`
        USEDPROC_CONST=`cat ${DGEMMETIMEDIR}/${N} | grep "ConstantPM->" | awk '{print $9}'`
				
        if [ $HETTIME ]; then
        	echo ${N} ${TIMETH} ${HETTIME} ${USEDPROC_HET} ${SMOOTHTIME} ${USEDPROC_SMOOTH} ${CONSTTIME} \
        				${USEDPROC_CONST}
      	else
      		echo " "
      	fi
done

#################################################################

exit 0

#################################################################

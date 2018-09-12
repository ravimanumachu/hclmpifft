#!/bin/bash

#################################################################

BLASBASEDIR=Results/
DGEMMETIMEDIR=${BLASBASEDIR}/parts
mkdir -p ${DGEMMETIMEDIR}

#################################################################

for (( COUNT=1; COUNT<=4000; COUNT+=1 ));
do
        N=`expr $COUNT \* 64`

				#echo ./workdist $N hetro 0 /home/hamidreza/hpopt/time_functions/set4/time

        ./workdist $N het 0 /home/hamidreza/hpopt/functions/set4/time/ > ${DGEMMETIMEDIR}/${N} 2>&1
        
        WTIME=`cat ${DGEMMETIMEDIR}/${N} | grep "Time Threshold=" | awk '{print $3}'`
        SOLUTIONTIME=`cat ${DGEMMETIMEDIR}/${N} | grep "ExecTime=" | awk '{print $3}'`
        USEDPROC=`cat ${DGEMMETIMEDIR}/${N} | grep "UsedProc=" | awk '{print $9}'`
        PARTITIONTIME=`cat ${DGEMMETIMEDIR}/${N} | grep "Total Execution Time" | awk '{print $10}'`
        NCALL=`cat ${DGEMMETIMEDIR}/${N} | grep "numCall" | awk '{print $3}'`
        MAXPOINT=`cat ${DGEMMETIMEDIR}/${N} | grep "Max Num Points=" | awk '{print $4}'`
        
        if [ $SOLUTIONTIME ]; then
        	echo ${N} ${WTIME} ${SOLUTIONTIME} ${USEDPROC} ${NCALL} ${PARTITIONTIME} ${MAXPOINT}
      	else
      		echo " "
      	fi
done

#################################################################

exit 0

#################################################################

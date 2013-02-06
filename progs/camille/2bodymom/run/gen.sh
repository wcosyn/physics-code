#!/bin/bash


n=3
arguments[0]="/user/home/gent/vsc407/vsc40761/Code/progs/camille/2bodymom/2bodymom_glauber" # the executable name
arguments[2]=$n # an extra parameter



# QSUB SETTINGS -----------------------------------------
walltime="walltime=00:01:00"
resource="nodes=1:ppn=1"
name="glauber_test_"
que="-q debug" # all extra qsub params go here together with option setters, ! carefull no space @ start or end !
# -------------------------------------------------------


echo "jobid: ${PBS_JOBID/.*/}" ## get the job id but strip everything after first .
echo "jobname: ${PBS_JOBNAME}"


workdir=`pwd` 		# location of the node.sh script
node_script=node.sh 	# name of node script
echo "" > runs.txt 	# clear the runs.txt file

for ((i=0; i<$n; i++ )) # C-style loop
do
	arguments[1]=$i
	echo "cd ${workdir}; bash node_script.sh ${arguments[@]}  | qsub ${que} -l ${resource} -l ${walltime} -N ${name}${i} " ##echo to terminal
	echo "cd ${workdir}; bash node_script.sh ${arguments[@]}  | qsub ${que} -l ${resource} -l ${walltime} -N ${name}${i} " >> "runs.txt" # echo to file
	echo "cd ${workdir}; bash node_script.sh ${arguments[@]}" | qsub ${que} -l "${resource}" -l "${walltime}" -N "${name}${i}" # really execute the thing
done

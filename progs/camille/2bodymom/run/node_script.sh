#!/bin/bash

function createtmpdir
{
	if [ $(hostname) = trillian ]
	then
		echo "Host trillian recognized"
		tmpdir="/tmp/`whoami`_`hostname`_${jobid}"
		echo "Creating dir ${tmpdir} "
		mkdir -p ${tmpdir}
	elif [[ $(hostname) = *.gengar.* ]] || [[ $(hostname) = *.gastly.* ]] || [[ $(hostname) = *.haunter.* ]] # double brackets to ensure * wildcard expansion
	then
		echo "Host gastly|haunter|gengar recognized"
		tmpdir="${VSC_SCRATCH_NODE}/${USER}_${jobid}_data"
		echo "Creating dir ${tmpdir}"
		mkdir -p ${tmpdir}
	else
		echo "Unrecognized host"
		exit
	fi
}

function run
{
	echo "$args > ${tmpdir}/stdout.log 2> ${tmpdir}/stderr.log" ## echo for log files
	$args > ${tmpdir}"/stdout.log" 2> ${tmpdir}"/stderr.log" # launch program, redirect output to logfiles
}

function copytempdir
{
	echo "contents of tmpdir ${tmpdir} : `ls -l ${tmpdir}`"
	echo "copying temp files to ${logdir}"
	mkdir -p ${logdir} # make logdir if not already existst
	for file in ${tmpdir}/*
	do
		mv $file "${logdir}/${PBS_JOBNAME}_${file/*\//}" # remove leading path from filenames -> /*\//
	done
}

# ================================================= ACTUAL SCRIPT ==========================================
# ==========================================================================================================
# ==========================================================================================================


jobid=${PBS_JOBID/.*/} ## get the job id but strip everything after first .
logdir="log"		# the dir of the log files; set this to your liking, (if relative -> rel to loc of this script)
args=$@			# received arguments are in $@ in this scope
argc=$#			# number of received arguments is in $#
# ======================================= some echoes to the log files =====================================

echo ""
echo "node script ${jobid} with job name ${PBS_JOBNAME} reporting with: $argc arguments received: $args"
echo "my current directory is `pwd`"
echo "log output directory is ${logdir}"
echo "hostname is: $(hostname)"
# ======================================= now finally do something =========================================

createtmpdir		# tmpdir is a now variable containing the path
run 			# run the program, program name and all variables are contained in $@ (passed via command line)
copytempdir		# now copy the contents of the tmpdir to our logdir
rm -r ${tmpdir}		# remove the tmpdir



# ==========================================================================================================
# ==========================================================================================================

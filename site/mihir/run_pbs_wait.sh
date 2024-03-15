#!/bin/bash
SELF=$(realpath $0)
PKGHOME=${SELF%/site*}
JOBSDIR="${PKGHOME}/jobs"

TYPE=${1##*.}
if [[ ${TYPE} == "sh" ]]; then
	export PYLAUNCH=$1
	shift
else
	export PYLAUNCH="${JOBSDIR}/py2launch.sh"
fi
echo ${PYLAUNCH}

export PYSCRIPT=$1
shift
export ARGS=$@

wait_flag=1
QUEUE_MAMU="NCMRWF1"
LOGDIR="/scratch/${USER}/logs"
PBSDIR="/scratch/${USER}/jobs"
PYFILE=${PYSCRIPT##*/}
TASKNAM=${PYFILE%.py}
PBSFILE="${PBSDIR}/run_task_${TASKNAM}.pbs"
RUNTIME=$(date +%Y%m%d_%H%M)

if [ ! -d ${LOGDIR} ]; then mkdir -p ${LOGDIR}; fi
if [ ! -d ${PBSDIR} ]; then mkdir -p ${PBSDIR}; fi

rm -f ${PBSFILE}
cat >> ${PBSFILE} << EOF
#!/bin/bash
#PBS -q ${QUEUE_MAMU}
#PBS -N ${TASKNAM}
#PBS -l select=1:ncpus=1:vntype=cray_compute
#PBS -o ${LOGDIR}/${TASKNAM}_${RUNTIME}.out
#PBS -e ${LOGDIR}/${TASKNAM}_${RUNTIME}.err
export LOGDIR=${LOGDIR}
export PBSDIR=${PBSDIR}
export TASKNAM=${TASKNAM}
export PKGHOME=${PKGHOME}


PYSCRIPT=${PYSCRIPT}
PYLAUNCH=${PYLAUNCH}
ARGS=${ARGS}

aprun -n 1 -N 1 \${PYLAUNCH} \${PYSCRIPT} \${ARGS}
EOF

cat $PBSFILE
jobid=$(qsub ${PBSFILE})

echo ${jobid}
	cnt=${wait_flag}
	while [[ ${cnt} -ne 0 ]]; do
		cnt=$(qstat -w ${jobid}|wc -l)
		qstat -w ${jobid}
		sleep 10
	done
tail -10 ${LOGDIR}/${TASKNAM}_${RUNTIME}.err
echo ${LOGDIR}
ls -lrt ${LOGDIR}/${TASKNAM}_${RUNTIME}.*

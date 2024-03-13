#!/bin/bash

export BASHJOB=$1
shift
export ARGS=$@

wait_flag=0
QUEUE_MAMU="NCMRWF1"
LOGDIR="/scratch/${USER}/logs"
PBSDIR="/scratch/${USER}/jobs"
JOBFILE=${BASHJOB##*/}
TASKNAM=${JOBFILE%.*}
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


BASHJOB=${BASHJOB}
ARGS=${ARGS}

\${BASHJOB} \${ARGS}
EOF

jobid=$(qsub ${PBSFILE})

echo ${jobid}
	cnt=${wait_flag}
	while [[ ${cnt} -ne 0 ]]; do
		cnt=$(qstat -w -u ${USER}|grep ${jobid}|wc -l)
	done

qstat -w -u ${USER}
echo ${LOGDIR}
ls -lrt ${LOGDIR}

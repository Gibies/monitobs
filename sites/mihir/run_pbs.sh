#!/bin/bash

export PYLAUNCH=$1
shift
export PYSCRIPT=$1
shift
export ARGS=$@

QUEUE_MAMU="serial1"
LOGDIR="/scratch/${USER}/logs"
PBSDIR="/scratch/${USER}/jobs"
PBSFILE="${PBSDIR}/run_task_vimt.pbs"
TASKNAM="vimt"
RUNTIME=$(date +%Y%m%d_%H%M)

if [ ! -d ${LOGDIR} ]; then mkdir -p ${LOGDIR}; fi
if [ ! -d ${PBSDIR} ]; then mkdir -p ${PBSDIR}; fi

rm -f ${PBSFILE}
cat >> ${PBSFILE} << EOF
#!/bin/bash
#PBS -q ${QUEUE_MAMU}
#PBS -N ${TASKNAM}
#PBS -l select=1:ncpus=1:vntype=cray_mamu
#PBS -o ${LOGDIR}/${TASKNAM}_${RUNTIME}.out
#PBS -e ${LOGDIR}/${TASKNAM}_${RUNTIME}.err
export LOGDIR=${LOGDIR}
export PBSDIR=${PBSDIR}
export TASKNAM=${TASKNAM}


PYSCRIPT=${PYSCRIPT}
PYLAUNCH=${PYLAUNCH}
ARGS=${ARGS}

\${PYLAUNCH} \${PYSCRIPT} \${ARGS}
EOF

jobid=$(qsub ${PBSFILE})

echo ${jobid}

qstat -w -u ${USER}|grep ${jobid}

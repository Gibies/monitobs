#!/bin/bash
SELF=$(realpath $0)
SITEBIN="${SELF%/*}"
SITEFLDR=${SITEBIN%/bin*}
SITE=${SITEFLDR##*/}
PKGHOME=${SELF%/site*}
PKGNAM=${PKGHOME##*/}
JOBSDIR="${PKGHOME}/jobs"

export JOBSCRIPT=$1
export HOSTPKG=${JOBSCRIPT%/jobs/*}
TYPE=${JOBSCRIPT##*.}
if [[ ${TYPE} == "py" ]]; then
	export JOBCMD="python ${JOBSCRIPT}"
	shift
else
	export JOBCMD=${JOBSCRIPT}
	shift
fi
export ARGS=$@

wait_flag=1
QUEUE_MAMU="NCMRWF1"
NODE_CNT=1
LOGDIR="/scratch/${USER}/logs"
PBSDIR="/scratch/${USER}/jobs"
PYFILE=${JOBSCRIPT##*/}
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
#PBS -l select=${NODE_CNT}:ncpus=1:vntype=cray_compute
#PBS -l place=scatter
#PBS -o ${LOGDIR}/${TASKNAM}_${RUNTIME}.out
#PBS -e ${LOGDIR}/${TASKNAM}_${RUNTIME}.err
export LOGDIR=${LOGDIR}
export PBSDIR=${PBSDIR}
export TASKNAM=${TASKNAM}
export PKGHOME=${PKGHOME}
export HOSTPKG=${HOSTPKG}
export NODE_CNT=${NODE_CNT}
export JOBSCRIPT=${JOBSCRIPT}
export JOBCMD="${JOBCMD}"
export ARGS=${ARGS}

export MODULEPATH="${PKGHOME}/site/${SITE}/modules:\${MODULEPATH}"
echo \${MODULEPATH}
module load public/py2env
module load ${USER}/${PKGNAM} 


export CUSLIB="\${HOSTPKG}/customlib"
if [ -d \${CUSLIB} ]; then 
export PYTHONPATH=\${CUSLIB}:\${PYTHONPATH}
fi

module list

which python

cd /scratch/${USER}
echo \${CUSLIB}

set -x
aprun -n 1 -N 1 \${JOBCMD} \${ARGS}
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
ln -sf ${LOGDIR}/${TASKNAM}_${RUNTIME}.out ${LOGDIR}/${TASKNAM}.out
ln -sf ${LOGDIR}/${TASKNAM}_${RUNTIME}.err ${LOGDIR}/${TASKNAM}.err
qstat -fx ${jobid}
cat ${LOGDIR}/${TASKNAM}.out
tail -10 ${LOGDIR}/${TASKNAM}.err
ls -lrt ${LOGDIR}/${TASKNAM}.*
echo ${LOGDIR}

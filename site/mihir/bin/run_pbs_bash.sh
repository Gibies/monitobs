#!/bin/bash
SELF=$(realpath $0)
SITEBIN="${SELF%/*}"
SITEFLDR=${SITEBIN%/bin*}
SITE=${SITEFLDR##*/}
PKGHOME=${SELF%/site*}
PKGNAM=${PKGHOME##*/}
JOBSDIR="${PKGHOME}/jobs"

TYPE=${1##*.}
if [[ ${TYPE} == "sh" ]]; then
	export BASHJOB=$1
	shift
	export HOSTPKG=${BASHJOB%/jobs/*}
else
	export BASHJOB="python $1"
	shift
fi
export ARGS=$@

wait_flag=1
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
export PKGHOME=${PKGHOME}
export HOSTPKG=${HOSTPKG}
export MODULEPATH="${PKGHOME}/site/${SITE}/modules:\${MODULEPATH}"
module load public/py2env
module load ${USER}/${PKGNAM} 

export CUSLIB="\${PYSCRIPT%/jobs/*}/customlib"
if [ ! -d \${CUSLIB} ]; then export CUSLIB="\${PKGHOME}/customlib"; fi
export PYTHONPATH=\${CUSLIB}:\${PYTHONPATH}

echo \${MODULEPATH}
module list

which python

BASHJOB=${BASHJOB}
ARGS=${ARGS}
aprun \${BASHJOB} \${ARGS}
EOF

cat ${PBSFILE}
jobid=$(qsub ${PBSFILE})

echo ${jobid}
	cnt=${wait_flag}
	while [[ ${cnt} -ne 0 ]]; do
		cnt=$(qstat -w ${jobid}|wc -l)
		qstat -w ${jobid}
		sleep 10
	done
echo ${LOGDIR}
ls -lrt ${LOGDIR}/${TASKNAM}_${RUNTIME}.*
tail -f ${LOGDIR}/${TASKNAM}_${RUNTIME}.err

#!/bin/bash
SELF=$(realpath ${0})
TASKNAME=${SELF##*/}
export GITROOT=${SELF%/site*}
export JOBSDIR=${GITROOT}/jobs
echo ${GITROOT}
cd ${GITROOT}
if [ -z ${HOST} ] ; then
${JOBSDIR}/${TASKNAME}
else
if [ ${HOST:0:6} == "elogin" ]; then
ssh -x utility01 ${JOBSDIR}/${TASKNAME}
else
${JOBSDIR}/${TASKNAME}
fi
fi

#!/bin/bash
SELF=$(realpath ${0})
HOST=$(hostname)
export JOBSDIR=${SELF%/*}
export GITROOT=${SELF%/jobs/*}

if [ ${HOST:0:6} == "elogin" ]; then
ssh -x utility01 'cd '${GITROOT}'; git pull ; git status'
else
cd ${GITROOT}; git pull ; git status
fi

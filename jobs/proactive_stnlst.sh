#!/bin/bash
module load gnu/pythonmodules/2.7.9/gnu/pandas/0.18.1
module load gnu/matplotlib/2.2.2

set -x
SELF=$(realpath $0)
CURRDIR=${SELF%/*}

PYSCRYPT="${CURRDIR}/read_fsoi.py"

TODATE=$(date +%Y%m%d)

export SUITE_ID=${SUITE_ID:-"Test"}
export ROSE_SUITE_DIR=${ROSE_SUITE_DIR:-"/scratch/${USER}/cylc-run/${SUITE_ID}"}
export CYLC_TASK_NAME=${CYLC_TASK_NAME:-"glu_proactive_stnlst"}
export VarNml=${VarNml:-"${ROSE_SUITE_DIR}/modules/monitobs/nml/varobs_nml"}
export StnLstDir=${StnLstDir:-"${ROSE_SUITE_DIR}/share/data/etc/stationlists/atmos"}
export CYLC_TASK_CYCLE_POINT=${CYLC_TASK_CYCLE_POINT:-"${TODATE}T0000Z"}
export INPUTDIR=${INPUTDIR:-"${ROSE_SUITE_DIR}/share/data/fsoimpact"}
export OUTDIR=${OUTDIR:-"${ROSE_SUITE_DIR}/share/cycle/${CYLC_TASK_CYCLE_POINT}/stationlists"}
export WRKDIR=${WRKDIR:-"${ROSE_SUITE_DIR}/work/${CYLC_TASK_CYCLE_POINT}/${CYLC_TASK_NAME}"}

python ${PYSCRYPT}

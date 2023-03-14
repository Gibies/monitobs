#!/bin/bash
SELF=$(realpath ${0})
export GITROOT=${SELF%/*}

ssh -x utility01 'cd '${GITROOT}'; git pull ; git status'

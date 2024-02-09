#!/bin/bash
SELF=$(realpath ${0})
echo ${SELF}
export JOBSDIR=${SELF%/*}
export GITROOT=${SELF%/jobs/*}
export HOST=$(hostname)
cd ${GITROOT}
UTFList=$(git ls-files --others --exclude-standard)
for utf in ${UTFList}; do
echo ${utf}
done
if [ ${HOST:0:6} == "elogin" ]; then
ssh -x utility01 'cd '${GITROOT}'; git add '${UTFList}' ; git status'
else
cd ${GITROOT}; git add ${UTFList} ; git status
fi
exit

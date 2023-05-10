#!/bin/bash
SELF=$(realpath ${0})
export GITROOT=${SELF%/*}

cd ${GITROOT}
UTFList=$(git ls-files --others --exclude-standard)

for utf in ${UTFList}; do
echo ${utf}
done
ssh -x utility01 'cd '${GITROOT}'; git add '${UTFList}' ; git status'

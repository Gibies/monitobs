#!/bin/bash
SELF=$(realpath ${0})
export GITROOT=${SELF%/*}

cd ${GITROOT}
UTFList=$(git ls-files --others --exclude-standard)

for utf in ${UTFList}; do
ssh -x utility01 'cd '${GITROOT}'; git add '${utf}' ; git status'
done

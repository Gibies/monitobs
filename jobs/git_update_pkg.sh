#!/bin/bash
set -x
SELF=$(realpath ${0})
export GITROOT=${SELF%/jobs*}
echo ${GITROOT}
cd ${GITROOT}
UTFList=$(git ls-files --others --exclude-standard)
####################################################
if [ ! -z "${UTFList}" ]; then
echo "_______________________"
echo "New untracked files are"
for utf in ${UTFList}; do
echo ${utf}
done
echo "_______________________"
fi
####################################################
cd ${GITROOT}; 
git pull
git add ${UTFList} 
git commit -a -m "Add files $(echo ${UTFList}| tr "
" " ")"
git push
git status



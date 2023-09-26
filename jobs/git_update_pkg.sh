#!/bin/bash
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
cd ${GITROOT}; 
git pull
git add ${UTFList} 
git commit -a -m "Add files $(echo ${UTFList}| tr "
" " ")"
git push
git status
####################################################
else
git pull
git commit -a -m "Minor change"
git push
git status
fi



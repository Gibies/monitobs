#!/bin/bash
set -x
SELF=$(realpath $0)
JOBSDIR=${SELF%/*}
HOMEDIR=${SELF%/jobs/*}
ROOTDIR=${SELF%/modules/*}

PKG_NAME=${PKG_NAME:-${HOMEDIR##*/}}

export PREFIX_PATH=${PREFIX_PATH:-$(pwd)}

while test $# -gt 0; do
echo $1 
        case "$1" in
                -h|--help)      shift;helpdesk;;
		-p|--prefix)	shift;PREFIX_PATH=$1;shift;;
                -u|--user)      shift;USER=$1;shift;;
                *)              ARGS="${ARGS} $1"; shift;;
        esac
done

#############################################################################################
#############################################################################################
if [ ! -d ${PREFIX_PATH} ]; then mkdir -p ${PREFIX_PATH}; fi

cd ${PREFIX_PATH}
pwd

sharedkey_list=""
indvdlkey_list="jobs src pylib pydic nml nglib"
key_list="${sharedkey_list} ${indvdlkey_list}"
for key in ${key_list}; do
if [ ! -d ${PREFIX_PATH}/${key} ]; then mkdir -p ${PREFIX_PATH}/${key} ; fi
if grep -q "${key}" <<< "${sharedkey_list}"; 
then
  echo "${PREFIX_PATH}/${key}"
  for file in $(find ${HOMEDIR}/${key}/*); do
  file_target=${PREFIX_PATH}/${key}/${file##*/${key}/}
  if [ ! -e ${file_target} ]; then
  echo ${file} ${file_target}
  cp -r ${file} ${file_target}
  fi
  done
fi
if grep -q "${key}" <<< "${indvdlkey_list}"; 
then
  echo "${PREFIX_PATH}/${key}"
  for file in $(find ${HOMEDIR}/${key}/*); do
  file_target=${PREFIX_PATH}/${key}/${file##*/${key}/}
  echo ${file} ${file_target}
  cp -r ${file} ${file_target}
  done
fi
done
#############################################################################################
#############################################################################################
set -x
NCOS_ROOT=${PREFIX_PATH%/SOURCE_DIR*}
SUITE_ID=${NCOS_ROOT##*/}
ROSE_APP_LIB=${NCOS_ROOT}/roses/${SUITE_ID}_Hybrid/app
if [ ! -d ${ROSE_APP_LIB} ]; then mkdir -p ${ROSE_APP_LIB}; fi
for app_path in $(ls -d ${HOMEDIR}/app/*_${obsname}*/); do
app_path=${app_path%/}
appname=${app_path##*/app/}
echo ${appname}
if [ -d ${ROSE_APP_LIB}/${appname} ]; then 
cd ${ROSE_APP_LIB}
APP_ARCHIVE="${ROSE_APP_LIB}/${appname}_backup_$(date +%Y%m%d_%H%M).tar.gz"
tar -cvzf ${APP_ARCHIVE} ${appname}
echo "${APP_ARCHIVE} is the backup"
rm -rf ${ROSE_APP_LIB}/${appname}; 
fi
cp -r ${app_path} ${ROSE_APP_LIB}/
sed -i "s+^OPPS_HOME=.*+OPPS_HOME=${PREFIX_PATH}+" ${ROSE_APP_LIB}/${appname}/rose-app.conf
cat ${ROSE_APP_LIB}/${appname}/rose-app.conf
done

#############################################################################################
#############################################################################################
set +x
export BACKUP="${PREFIX_PATH}/jobs/packup_${PKG_NAME}.sh"
if [ -f ${BACKUP} ]; then 
chmod 750 ${BACKUP}
${BACKUP} ; 
else
echo "${BACKUP} file not found"
fi
#############################################################################################

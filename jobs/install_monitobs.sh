#!/bin/bash
set -x
SELF=$(realpath $0)
JOBSDIR=${SELF%/*}
HOMEDIR=${SELF%/jobs/*}
ROOTDIR=${SELF%/modules/*}

PKG_NAME=${PKG_NAME:-${HOMEDIR##*/}}


while test $# -gt 0; do
echo $1 
        case "$1" in
                -h|--help)      shift;helpdesk;;
		-m|--modpath)	shift;PKGSYSMOD=$1;shift;;
		-p|--prefix)	shift;PREFIX_PATH=$1;shift;;
                -u|--user)      shift;USER=$1;shift;;
		-v|--version)	shift;VRSN=$1;shift;;
                *)              ARGS="${ARGS} $1"; shift;;
        esac
done

export VRSN=${VRSN:-$(date +%Y.%m.%d)}
export VRNT=${VRNT:-$(date +%H%M)}
export PREFIX_PATH=${PREFIX_PATH:-$(pwd)/${PKG_NAME}.${VRSN}}
if [ ! -d ${PREFIX_PATH} ]; then mkdir -p ${PREFIX_PATH}; fi
echo "${PKG_NAME} ${VRSN} ${VRNT}" > ${PREFIX_PATH}/version
#############################################################################################
#############################################################################################

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
set +x
export BACKUP="${PREFIX_PATH}/jobs/packup_${PKG_NAME}.sh"
if [ -f ${BACKUP} ]; then 
chmod 750 ${BACKUP}
${BACKUP} ; 
else
echo "${BACKUP} file not found"
fi
#############################################################################################
export PKGSYSMOD=${PKGSYSMOD:-${PREFIX_PATH}/modules}
if [ ! -d ${PKGSYSMOD} ]; then mkdir -p ${PKGSYSMOD}; fi
export modulefile=${PKGSYSMOD}/${PKG_NAME}.${VRSN}
rm -f ${modulefile}
cat >> ${modulefile} << EOF
#%Module1.0###
###
### modules modulefile
###
proc ModulesHelp {  } {
puts stderr "${modulefile}"
}
module-whatis   "${PKG_NAME}.${VRSN}"
set     TOP     ${PREFIX_PATH}
setenv	${PKG_NAME} \$TOP
setenv	${PKG_NAME}_lib	\$TOP/pylib
setenv	${PKG_NAME}_dic	\$TOP/pydic
setenv	${PKG_NAME}_nml	\$TOP/nml
prepend-path    PATH    \$TOP/jobs
prepend-path    PYTHONPATH    \$TOP/pylib:\$TOP/pydic:\$TOP/nml

if { [ module-info mode load ] } {
puts stdout ""
}
EOF

set -x
export MODULEPATH=$MODULEPATH:${PKGSYSMOD}
module show ${modulefile}
#############################################################################################
export sample="${PREFIX_PATH}/jobs/sample.py"
rm -f ${sample}
cat >> ${sample} << EOF
import os,sys
sys.path.append(os.environ.get('monitobs_lib'))
import modulib

print("Test is done")
EOF
#############################################################################################
${PREFIX_PATH}/jobs/py2launch.sh ${sample}

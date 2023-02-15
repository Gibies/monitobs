#!/bin/bash
set -x
SELF=$(realpath $0)
CURRDIR=${SELF%/*}
HOMEDIR=${SELF%/jobs/*}
SUITE_DIR=${SELF%/modules/*}
PKG_NAME=${PKG_NAME:-${HOMEDIR##*/}}
WRKDIR=${WRKDIR:-"/scratch/${USER}/packages/${PKG_NAME}"}
TARBALL=${TARBALL:-${1:-"${WRKDIR}.tar.gz"}}

if [ -d ${WRKDIR} ]; then rm -rf ${WRKDIR} ; fi
mkdir -p ${WRKDIR} ;
mkdir -p ${WRKDIR}/modules
mkdir -p ${WRKDIR}/app

cp -rL ${HOMEDIR} ${WRKDIR}/modules
for appname in $(cat ${HOMEDIR}/nml/applist); do
cp -rL ${SUITE_DIR}/app/${appname} ${WRKDIR}/app
done

cd ${WRKDIR}
tar -zcvf ${TARBALL} modules app
rm -rf ${WRKDIR}
echo ${TARBALL}

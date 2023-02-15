#!/bin/bash
runstart=`date +%s`
#################################################################################
LIBDIR=/home/apps/SiteSoftwares/intel/NWPROD-LIB/bufr-intel/11.2.0
module switch PrgEnv-cray/6.0.4 PrgEnv-intel/6.0.4
module load intel/eccodes/2.9.2_1
module load gnu/pythonmodules/2.7.9/gnu/pandas/0.18.1
module load gnu/utility/pynio/1.5.5
module load gnu/pyngl/1.6.1
module load gnu/pythonmodules/2.7.9/gnu/netcdf4/1.5.3-netcdf4.6.0hdf51.10.0
module load gnu/matplotlib/2.2.2
module load gnu/pythonmodules/2.7.9/gnu/basemap/1.1.0
module load gnu/pythonmodules/2.7.9/gnu/py-ncepbufr/1.1.1
module load gnu/python_eccodes/eccodes-2.13.1_utility
module load /home/${USER}/SoftwarePackages/modules/aapp-8.9
#################################################################################

SELF=$(realpath ${0})
export HOMEDIR=${SELF%/jobs/*}
obsname="mwri"
if [ ${USER} == "umprod" ]; then 
export QUEUE_MAMU="serial"
else
export QUEUE_MAMU="serial1"
fi
export SCRATCH="/scratch/${USER}/tmpFldr/obstore"
export BUFRTNK="/home/${USER}/data/bufr_tank"
#DataSorcKey="_gfsprod_gdas_dump"
#export GDASDIR="/home/gfsprod/data/gdasv14/gdas/prod"	#/gdas.20210225/gdas1.t06z.satwnd.tm00.bufr_d
#DataSorcKey="_umprod_gdas_dump"
#export GDASDIR="/home/${USER}/data/gfsprod/gdas/prod"	#/gdas.20210225/gdas1.t06z.satwnd.tm00.bufr_d
DataSorcKey="ecbufr"
#obstype="20700"
#elemcnt="19"
export JOBDIR=${HOMEDIR}/jobs
export SRPTDIR=${HOMEDIR}/scripts
export SRCDIR=${HOMEDIR}/src
export NMLDIR=${HOMEDIR}/nml
#export obsname="$(grep "obsname" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export gdasname="$(grep "gdasname" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export bufrmtyp="$(grep "bufrmtyp" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export obsubtyp="$(grep "obsubtyp" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export elemcnt="$(grep "elemcnt" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export OBSTORE="$(grep "OBSTORE" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export obsgroup="$(grep "obsgroup" ${NMLDIR}/keys_${obsname}.nml | awk '{print $2}' )"
export WRKDIR=${WRKDIR:-${SCRATCH}/work/${obsname}}
export LOGDIR=${LOGDIR:-${SCRATCH}/log/${obsname}}
export EXEDIR=${EXEDIR:-${WRKDIR}/exec}
export WRKHDRDIR=${WRKHDRDIR:-${HOMEDIR}/hdr}
export HDR_FILE="${obsname}.hdr"     ### Name of header file
#################################################################################
while test $# -gt 0; do
echo $1 $2
        case "$1" in
		-c|--compile)	${JOBDIR}/compile_opps_${obsname}.sh ; shift;;
                -d|--date)      shift; DATE_TIME=$(date -d "$1" +"%Y%m%d %H:%M") ; shift;;
                -h|--help)      helpdesk; shift;;
		-t|--windowth)	shift;export WNDWHALF=$1; shift;;
                *)              ARGS="${ARGS} $1"; shift;;
        esac
done
#################################################################################
set -x
if [ -z ${DATE_TIME} ]; then 
PDY=${PDY:-$(date +%Y%m%d)}
CYC=${CYC:-"00"}
DATE_TIME=$(date -d "${PDY} ${CYC}" +"%Y%m%d %H:%M"); 
else
PDY=${PDY:-$(date --date="${DATE_TIME}" +%Y%m%d)}
CYC=${CYC:-$(date --date="${DATE_TIME}" +%H)}
fi
if [ -z ${WNDWHALF} ]; then WNDWHALF=3; fi
RDATE=$(date --date=${PDY} +%d%m%Y)

MINUTES_WINDOW_HALF=$(echo "(${WNDWHALF} * 60) - 1" | bc -l)
START_TIME=$(date -d "${DATE_TIME} ${MINUTES_WINDOW_HALF} minutes ago" +"%Y%m%d %H:%M")
END_TIME=$(date -d "${DATE_TIME} ${MINUTES_WINDOW_HALF} minutes" +"%Y%m%d %H:%M")
echo ${START_TIME} ${END_TIME}

curtime=$(date -d "${START_TIME}" +"%Y%m%d%H%M")
endtime=$(date -d "${END_TIME}" +"%Y%m%d%H%M")
echo ${curtime} ${endtime}
#################################################################################
export LOGOUTFILE=${LOGDIR}/${obsname}.log
if [ ! -d ${LOGDIR} ]; then mkdir -p ${LOGDIR}; fi
#exec &> >(tee -a ${LOGOUTFILE} >&1 )
#################################################################################
export WRKBUFRDIR=${WRKBUFRDIR:-${WRKDIR}/${PDY}/${CYC}/workdir_bufr}
export WRK_OBSTORE=${WRK_OBSTORE:-${WRKDIR}/${PDY}/${CYC}/workdir_obstore}
export OUTDIR=${OUTDIR:-${SCRATCH}/output/${obsname}${DataSorcKey}_${PDY}_${CYC}}
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; fi
if [ ! -d ${WRKDIR} ]; then mkdir -p ${WRKDIR} ; fi
#################################################################################
#################################################################################
#Section 4
##################Listing of the BUFR file#######################################
export BUFRSORCDIR="/home/gfsprod/data/idata/eumet"
export BUFRFLDRSTR="DBNet/mwri"
export BUFRFPREFIX="W_XX-EUMETSAT-*_C_EUMS_*_mwri_*"
export BUFRPOSTFIX="*_l1b_bufr.bin"
if [ -d ${WRKBUFRDIR} ]; then rm -rf ${WRKBUFRDIR}; fi
mkdir -p ${WRKBUFRDIR}
if [ -d ${WRK_OBSTORE} ]; then rm -rf ${WRK_OBSTORE}; fi
mkdir -p ${WRK_OBSTORE}

while [[ ${curtime} -le ${endtime} ]]; do
YYYY=${curtime:0:4}
MM=${curtime:4:2}
DD=${curtime:6:2}
HH=${curtime:8:2}
echo ${YYYY}_${MM}_${DD}_${HH}
YYYYMMDD=${YYYY}${MM}${DD}
YYYYMMDDHH=${YYYYMMDD}${HH}

curtime=$(date -d "${YYYY}${MM}${DD} ${HH} 1 hour" +"%Y%m%d%H%M")
cp -r ${BUFRSORCDIR}/${YYYYMMDD}/${BUFRFLDRSTR}/${BUFRFPREFIX}${YYYYMMDD}_${HH}*${BUFRPOSTFIX} ${WRKBUFRDIR}

done
ls -lrt ${WRKBUFRDIR}
#################################################################################
#################################################################################
#Section 5	AAPP is utilised for derived fields
#################################################################################
cd ${WRKBUFRDIR}
for file in $(ls ${WRKBUFRDIR}); do

echo $file
eccodes_decodebufr_1c -i ${file} MWRI

  export ORIGINATING_CENTRE=74
  export MASTER_TABLE=15
  export LOCAL_TABLE=1
  export MESSAGE_SUBTYPE=61
  export MWRI_THIN=-2
aapp_encodebufr_1c -i ${file%.*}.l1c MWRI1D

done
ls -lrt ${WRKBUFRDIR}
#################################################################################
PYSCRYPT="${SRPTDIR}/bufr_to_obstore_eccode.py"

if [ ! -d ${SRPTDIR} ]; then mkdir -p ${SRPTDIR} ; fi
rm -rf ${PYSCRYPT}
cat >> ${PYSCRYPT} << EOF

#################################################################################
from __future__ import print_function
import traceback
from eccodes import *

import os,sys

USER=os.environ.get('USER',"myhome")
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
print(OBSLIB)
import obsmod

if len(sys.argv) > 1:
    PDY = sys.argv[1]
else:
    PDY = obsmod.today()

if len(sys.argv) > 2:
    CYC = sys.argv[2]
else:
    CYC = "00"

obsname=os.environ.get('obsname', "${obsname}")
cntmax=int(os.environ.get('OBS_CNT_MAX', 500000))
btchcnt=int(os.environ.get('BTCH_CNT_MAX',5))
TDATE=os.environ.get('PDY',PDY)
print(TDATE)
Tnode=obsmod.pydate(TDATE)
datacfldrname=obsmod.cylcdate(Tnode)
outpath=os.environ.get('WRK_OBSTORE',"/scratch/"+USER+"/obstore/work/"+obsname+"/"+PDY+"/"+CYC+"/workdir_obstore")
inpath=os.environ.get('BUFRDIR',"")
slctstr=os.environ.get('BUFRFILESTR',"")
nmlfile=os.environ.get('KEYNMLFILE',OBSNML+"/keys_"+obsname+".nml")
eleindxmaptbl=os.environ.get('ELEINDXMAPTBL',OBSNML+"/aapp_fieldname.nml")

data=obsmod.ecbufr_decode_files(inpath,Tnode,slctstr,nmlfile,eleindxmaptbl=eleindxmaptbl)
obsmod.ascii_file_write(data,outfile=outpath+"/data_"+obsname+".txt",option=1)
outfile=obsmod.obstore_write(data,nmlfile,outpath,btchcnt=btchcnt,cntmax=cntmax,DT=Tnode,diagflag=0)

#################################################################################

EOF

set -x
export BUFRDIR=${WRKBUFRDIR}
export BUFRFILESTR="${BUFRFPREFIX}*${BUFRPOSTFIX%.*}.bufr"
export ELEINDXMAPTBL="${HOMEDIR}/nml/aapp_fieldname.nml"
export KEYNMLFILE="${NMLDIR}/keys_${obsname}.nml"
export BTCH_CNT_MAX=$(ls ${BUFRDIR}/${BUFRFILESTR} | wc -l)
if [[ ${BTCH_CNT_MAX} -gt 10 ]] ; then BTCH_CNT_MAX=10 ; fi
export OBS_CNT_MAX=200000
python ${PYSCRYPT} ${PDY} ${CYC}

wait
#################################################################################
#################################################################################
#Section 1
#################Compiler is called if Executables are not available#############
if [ ! -f ${EXEDIR}/obstore_${obsname}.exe ]; then
${JOBDIR}/compile_${obsname}.sh
fi
#################################################################################
#Section 2
#################Header will be created only if not available####################
if [ ! -d ${WRKHDRDIR} ]; then mkdir -p ${WRKHDRDIR}; fi
cd ${WRKHDRDIR}
pwd
if [ ! -f ${WRKHDRDIR}/${HDR_FILE} ]; then
${EXEDIR}/hdr_${obsname}.exe
echo "New Header file ${WRKHDRDIR}/${HDR_FILE} is created"
else
echo "Header file ${WRKHDRDIR}/${HDR_FILE} already exist"
fi
if [ ! -f ${WRKHDRDIR}/${HDR_FILE} ]; then
echo "Header file creation failed"
exit
fi

#################################################################################
#################################################################################
#Section 5
##################Obstore preperation############################################
cd ${WRK_OBSTORE}
pwd
export FIXHDR="${HOMEDIR}/fixhdr/fixhdr_608.hdr"
##########################################
#rm -f ${WRK_OBSTORE}/${OBSTORE}.obstore
##########################################
#cat ${WRK_OBSTORE}/${OBSTORE}.obstore.dat > ${WRK_OBSTORE}/${obsubtyp}.dat
#cp ${FIXHDR} ${WRK_OBSTORE}
#cp ${WRKHDRDIR}/${HDR_FILE} ${WRK_OBSTORE}
#ls ${WRK_OBSTORE}
#echo ${CYC}${RDATE} >ukmodate
#date --date="$PDY" +%j>>ukmodate
#date +"%Y%m%d%k%M%S" >>ukmodate
#date +"%j"   >>ukmodate
#cat ./ukmodate
#ls -al ukmodate
#set -x
#wc -l ${obsubtyp}.dat | awk '{print $1}' > cnt
#echo ${elemcnt}>> cnt
#${EXEDIR}/totobs_${obsname}.exe
#${EXEDIR}/obstore_${obsname}.exe
#################################################################################
### Section 6 : Copy obstore file to destination
#################################################################################
#WRKDIR_OUT_FILE="${WRK_OBSTORE}/${OBSTORE}.obstore"
#WRKDIR_OUT_FILE="${WRK_OBSTORE}/obstore.data"
#################################################################################
#################################################################################
WRKDIR_OUT_FILE="${WRK_OBSTORE}/${OBSTORE}.obstore"
FINAL_OUT_FILE="${OUTDIR}/${OBSTORE}.obstore"
if [[ -e ${WRKDIR_OUT_FILE} ]]; then
mv ${WRKDIR_OUT_FILE} ${FINAL_OUT_FILE}
echo ${FINAL_OUT_FILE} > ${WRK_OBSTORE}/outfile
echo "${OBSNAME} OBSTORE Completed Successfully..."
else
echo "${OBSNAME} OBSTORE Creation failed"
fi
#################################################################################
runend=`date +%s`
runtime=$((runend-runstart))
if [ $runtime -gt 60 ]; then
runminutes=$(echo "${runtime}/60" |bc -l)
echo "Runtime is ${runminutes} minutes"
else
echo "Runtime is ${runtime} seconds"
fi
exit
exit



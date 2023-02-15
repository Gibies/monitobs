#!/bin/sh
####################################################################################################################################
helpdesk()
{
                        echo "Options for "$0":"
			echo "___________________"
			echo "-d, --date "YYYYMMDD HH:MM"		specify the date and time"
			echo "-diag, --diag				diagnostics mode"
			echo "-dsp, --display				display the final pdf document"
                        echo "-h, --help               			show brief help"
                        echo "-t, --windowth				forecast window"
                        exit 0
}
####################################################################################################################################

####################################################################################################################################
options()
{
while test $# -gt 0; do
echo $1 $2
        case "$1" in
                -d|--date)      shift; DATE_TIME=$(date -d "$1 $2" +"%Y%m%d %H:%M") ; shift;shift;;
		-dig|--diag)	shift;export diag_flag="True";;
		-dsp|--display)	shift;export disp_flag="True";;
                -e|--edate)      shift; END_DATE_TIME=$(date -d "$1 $2" +"%Y%m%d %H:%M") ; shift;shift;;
                -h|--help)      helpdesk; shift;;
		-t|--windowth)	shift;export LEADTIME=$1; shift;;
                *)              ARGS="${ARGS} $1"; shift;;
        esac
done
}
####################################################################################################################################
SELF=$(realpath $0)
export HOMEDIR=${SELF%/jobs/*}
export JOBSDIR=${HOMEDIR}/jobs
export PYTHON_DIR="${HOMEDIR}/jobs"
export PKGNAME=${HOMEDIR##*/}

disp_flag="True"
TODATE=$(date +%Y%m%d)
options $(echo $@  | tr "=" " ")
echo ${DATE_TIME} ${END_DATE_TIME} ${LEADTIME}
####################################################################################################################################
if [[ -z ${LEADTIME} ]]; then LEADTIME=6 ; fi
if [[ -z ${DATE_TIME} ]]; then
DATE_TIME=$(date -d "1 day ago ${TODATE}" +"%Y%m%d %H:%M")
fi
if [[ -z ${END_DATE_TIME} ]]; then
if [[ -z ${NUMCYCLE} ]]; then
NUMCYCLE=1
fi
HOURFWRD=$(echo "${LEADTIME}*(${NUMCYCLE}-1)"| bc -l )
END_DATE_TIME=$(date -d "${HOURFWRD} hour ${DATE_TIME}" +"%Y%m%d %H:%M")
fi

CYLC_TIME_LIST=""

###########################################################################################
curtime=$(date -d "${DATE_TIME}" +"%Y%m%d%H%M")
endtime=$(date -d "${END_DATE_TIME}" +"%Y%m%d%H%M")
while [[ ${curtime} -le ${endtime} ]]; do
YYYY=${curtime:0:4}
MM=${curtime:4:2}
DD=${curtime:6:2}
HH=${curtime:8:2}
echo ${YYYY}_${MM}_${DD}_${HH}
YYYYMMDD=${YYYY}${MM}${DD}
CYLC_TIME=$(date -d "${YYYYMMDD} ${HH}" +%Y%m%dT%H00Z)
export CYLCDATE=$(echo ${CYLC_TIME} | cut -c1-8)
export CYLCHOUR=$(echo ${CYLC_TIME} | cut -c10-11)
###########################################################################################
CYLC_TIME_LIST="${CYLC_TIME_LIST} ${CYLC_TIME}"
####################################################################################################################################
curtime=$(date -d "${LEADTIME} hour ${YYYY}${MM}${DD} ${HH}" +"%Y%m%d%H%M")
done
###########################################################################################

set -x
INDATA_PATH="/home/umprod/cylc-run/op_trunk/share/cycle"
OUTDATA_PATH="/home/${USER}/data/fsoi/satwind"
if [[ ! -d ${OUTDATA_PATH} ]] ; then mkdir -p ${OUTDATA_PATH} ; fi

#CYLC_TIME_LIST="20220822T0000Z " #20211110T0600Z 20211110T1200Z 20211110T1800Z 20211111T0000Z 20211111T0600Z 20211111T1200Z 20211111T1800Z 20211112T0000Z"
subtype_list="23613 " #23632 23623 23621 23501 23502 23503 22501 22502 22505 22511 25203 25205 25207 25211 25202 25901 25902 25903"

which python
for CYLC_TIME in ${CYLC_TIME_LIST}; do

#for subtype in ${subtype_list}; do


echo "${subtype}_${CYLC_TIME}"
export subtype=${subtype}
export CYLC_TASK_CYCLE_POINT=${CYLC_TIME}
export INPUTDIR=${INDATA_PATH}
export OUTPUTDIR=${OUTDATA_PATH}
export INFILE="${INDATA_PATH}/${CYLC_TIME}/*.FSO"
export OUTFILE="${OUTDATA_PATH}/data_${subtype}_${CYLC_TIME}.fso"
${HOMEDIR}/jobs/pylaunch.sh ${PYTHON_DIR}/read_fsoi.py 
#/usr/bin/grep -ir " ${subtype} " ${INDATA_PATH}/${CYLC_TIME}/*.FSO > ${OUTDATA_PATH}/data_${subtype}_${CYLC_TIME}.fso

#done
done

ls -lrt ${OUTDATA_PATH}

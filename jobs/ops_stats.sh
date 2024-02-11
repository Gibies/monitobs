
CYLCTIME=${1:-"20240210T1200Z"}

INPATH="/home/umprod/PS43/LOG_DA_RUN"
OUTFILE="/scratch/${USER}/ops_count_info_${CYLCTIME}.txt"



grep "Total" ${INPATH}/${CYLCTIME}/glm_ops_process_background_*/NN/job.stats|awk '{$1=""; print $0}' > ${OUTFILE}

ls -lrt ${OUTFILE}

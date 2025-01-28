#!/bin/bash

# Watch for new sequencing runs to launch DRAGEN BCL-Convert
# TODO: Extract stats when FastqComplete.
# USAGE: Launch in crontab
#        bash dragen_bcl-convert_watcher.sh | tee -a ${LOGFILE}
#        bash /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_watcher.sh | tee -a /staging2/dragen/dragen_bcl-convert_watcher.log

# Scan BCL output dirs (BASEDIR) for new runs (FC) to demux
# Skip runs for LowPass (check if SampleSheet exists (LowPass))
# Demux outputs are written to ${WORKDIR}/dragen/${FC}
# File ${BASEDIR}/${FC}/CopyComplete.txt marks end of sequencing run
# File ${BASEDIR}/${FC}/Failed.txt marks that sequencing has failed
# File ${BASEDIR}/${FC}/FastqComplete.txt marks end of BCL-conversion

BASEDIR='/mnt/spxp-app02/staging/hiseq_raw'
WORKDIR='/mnt/spxp-app02/staging2/dragen'
LOGFILE="${WORKDIR}/dragen_bcl-convert_watcher.log"
LOGPREFIX="[bcl-watcher]"
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/A00977" "${BASEDIR}/LH00207R" "/mnt/vs_nas_chusj/SPXP_APP02_NFS/LH00336")

#printf "\n\n######\n%s %s %s\n######\n\n" $0 ${LOGPREFIX} $( date "+%F@%T" ) | tee -a ${LOGFILE}/

launch_run() {
    # Check if sequencing is finished (CopyComplete.txt) and that run
    # is not already being processed by another instance of this script
    # which creates output dir and Sample* files through `get_nanuq_files.py``
    local dir="$1"
    local fc="$2"
    parts=($(echo ${fc} | tr '_' '\n'))
    fc_short="${parts[1]}_${parts[2]}"
    if [[ -f "${dir}/${fc}/CopyComplete.txt" ]]; then
        echo "${LOGPREFIX} CopyComplete.txt file indicates that sequencing has finished"
        if [[ -d ${WORKDIR}/${fc} ]]; then
            echo "${LOGPREFIX} PASS: Run already processed in ${WORKDIR}/${fc}"
        else
            mkdir ${WORKDIR}/${fc}
            cd ${WORKDIR}/${fc}
            echo "${LOGPREFIX} Getting SampleSheet and other files from Nanuq..."
            # if [[ ! -f ${WORKDIR}/${FC} ]]; then
            # TODO: Skip get_nanuq_files if SampleSheet already present?
            # fi
            python /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py --run ${fc_short}
            if [[ -f "${WORKDIR}/${fc}/SampleSheet.csv" ]]; then
                echo "${LOGPREFIX} RUN: Launching BCL-convert with qsub..."
                qsub /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_launcher.sh ${fc}
            else
                echo "${LOGPREFIX} ERROR: SampleSheet.csv not found in ${WORKDIR}/${fc}" >&2
                exit 1
            fi
        fi
    else
        echo "${LOGPREFIX} PASS: Sequencing not finished. Waiting for CopyComplete.txt"
    fi
}

for dir in ${WATCHDIRS[@]}; do
    echo "${LOGPREFIX} Scanning ${dir}..."
    for FC in $( ls ${dir} ); do
        parts=($(echo ${FC} | tr '_' '\n'))
        if [[ ${#parts[@]} -eq 4 ]]; then
            # Check SampleSheet if run is LowPass.
            # If no SampleSheet is found, then run is not LowPass
            echo "---------------------------------"
            echo "${LOGPREFIX} WATCH: ${FC}"
            if [[ -f "${dir}/${FC}/SampleSheet.csv" ]]; then
                if grep -q "LowPass" "${dir}/${FC}/SampleSheet.csv"; then
                    echo "${LOGPREFIX} PASS: Found the word LowPass in SampleSheet"
                else
                    echo "${LOGPREFIX} WATCH: SampleSheet exists and not for LowPass."
                    launch_run ${dir} ${FC}
                fi
            elif [[ -f "${dir}/${FC}/LowPass*.csv" ]]; then
                echo "${LOGPREFIX} PASS: Found what looks like a LowPass SampleSheet."
            elif [[ -f "${dir}/${FC}/Failed.txt" ]] ||  [[ -f "${dir}/${FC}/failed.txt" ]]; then
                echo "${LOGPREFIX} PASS: Failed.txt marks a failed Run."
            elif [[ -f "${dir}/${FC}/FastqComplete.txt" ]]; then
                echo "${LOGPREFIX} PASS: FastqComplete.txt indicates that run has already been processed."
            else
                echo "${LOGPREFIX} WATCH: Could not find ${dir}/${FC}/SampleSheet.csv."
                launch_run ${dir} ${FC}
            fi
        fi
        # else: ignore because format of folder name doesn't look like a run
    done
    echo "================================="
done
exit 0

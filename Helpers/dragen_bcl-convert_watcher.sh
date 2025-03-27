#!/bin/bash

# Watch for new sequencing runs to launch DRAGEN BCL-Convert
# USAGE: Scheduled launch in /etc/cron.d/dragen_bcl-convert_watcher (see below)
#        bash dragen_bcl-convert_watcher.sh | tee -a ${LOGFILE}
#        bash /staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_watcher.sh | tee -a /staging2/dragen/dragen_bcl-convert_watcher.log

# /etc/cron.d/dragen_bcl-convert_watcher file content:
# # Run dragen_bcl-convert_watcher.sh and rotate its' log every Sunday
# */30 * * * * root /usr/bin/bash /staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_watcher.sh >> /mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/dragen_bcl-convert_watcher.log
# 59 23 * * Sun root path='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/dragen_bcl-convert_watcher'; flog="${path}.log"; dlog="${path}-$( date +"%Y%m%d" ).log"; mv ${flog} ${dlog}; gunzip ${flog}.tar.gz; tar -rf ${flog}.tar ${dlog}; gzip ${flog}.tar; rm -f ${dlog};

# Scan BCL output dirs (BASEDIR) for new runs (FC) to demux
# Skip runs for LowPass (check if SampleSheet exists (LowPass))
# Demux outputs are written to ${WORKDIR}/dragen/${FC}
# File ${BASEDIR}/${FC}/CopyComplete.txt marks end of sequencing run
# File ${BASEDIR}/${FC}/Failed.txt marks that sequencing has failed
# File ${BASEDIR}/${FC}/FastqComplete.txt marks end of BCL-conversion
#   FastqComplete.txt is copied to ${BASEDIR} by dragen_bcl-convert_launcher.sh

# DEPENDENCIES:
#   /staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_launcher.sh
#   /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py
# TODO: Extract stats when FastqComplete.

BASEDIR='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs'
HISEQ_R='/mnt/spxp-app02/staging/hiseq_raw'
WORKDIR='/mnt/spxp-app02/staging2/dragen'
LOGFILE="${WORKDIR}/dragen_bcl-convert_watcher.log"
LOGPREFIX="[bcl-watcher]"
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/LH00207R" "${HISEQ_R}/LH00336" "${HISEQ_R}/A00977" "${HISEQ_R}/LH00207R")

printf "\n\n######\n%s %s %s\n######\n\n" ${LOGPREFIX} $( date "+%F@%T" ) $0 #| tee -a ${LOGFILE}/

launch_run() {
    # Run dragen_bcl-convert_launcher.sh if not already being processed by 
    # another instance of this script which creates output dir ${WORKDIR}/${fc}
    # FastqComplete.txt is copied to ${BASEDIR} by dragen_bcl-convert_launcher.sh
    # and marks the run as done.
    local dir="$1"
    local fc="$2"
    parts=($(echo ${fc} | tr '_' '\n'))
    fc_short="${parts[1]}_${parts[2]}"
    if [[ -d "${WORKDIR}/${fc}/1.fastq" ]]; then
        echo "${LOGPREFIX} PASS: Demux appears to be in progress for ${WORKDIR}/${fc}"
    else
        mkdir ${WORKDIR}/${fc}
        cd ${WORKDIR}/${fc}
        echo "${LOGPREFIX} Getting SampleSheet and other files from Nanuq..."
        python /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py --run ${fc_short}
        if [[ -f "${WORKDIR}/${fc}/SampleSheet.csv" ]]; then
            echo "${LOGPREFIX} RUN: Launching BCL-convert with qsub..."
            qsub /staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_launcher.sh ${fc}
        else
            echo "${LOGPREFIX} ERROR: SampleSheet.csv not found in ${WORKDIR}/${fc}" >&2
        fi
    fi
}

for dir in ${WATCHDIRS[@]}; do
    echo "${LOGPREFIX} Scanning ${dir}..."
    for fc in $( ls ${dir} ); do
        parts=($(echo ${fc} | tr '_' '\n'))
        if [[ ${#parts[@]} -eq 4 ]]; then
            echo "---------------------------------"
            echo "${LOGPREFIX} ${fc}"
            if [[ -f "${dir}/${fc}/CopyComplete.txt" ]]; then
                echo "${LOGPREFIX} CopyComplete.txt indicates that sequencing has finished"
                # Check if bcl-convert needed (not previously demuxed, not failed, not LowPass)
                if [[ -f "${dir}/${fc}/FastqComplete.txt" ]]; then
                    echo "${LOGPREFIX} PASS: FastqComplete.txt indicates that run has already been processed."
                elif [[ -f "${dir}/${fc}/Failed.txt" ]] ||  [[ -f "${dir}/${fc}/failed.txt" ]]; then
                    echo "${LOGPREFIX} PASS: Failed.txt marks a failed Run."
                elif [[ -f "${dir}/${fc}/LowPass*.csv" ]]; then
                    echo "${LOGPREFIX} PASS: Found what looks like a LowPass SampleSheet."
                else
                    # Check SampleSheet if run is LowPass.
                    # If no SampleSheet is found, then run is not LowPass
                    if [[ -f "${dir}/${fc}/SampleSheet.csv" ]]; then
                        if grep -q "LowPass" "${dir}/${fc}/SampleSheet.csv"; then
                            echo "${LOGPREFIX} PASS: Found the word LowPass in SampleSheet"
                        elif grep -q "Cloud_Workflow," "${dir}/${fc}/SampleSheet.csv"; then
                            echo "${LOGPREFIX} PASS: SampleSheet indicates a Cloud_Workflow"
                        else
                            echo "${LOGPREFIX} SampleSheet exists and not for LowPass."
                            launch_run ${dir} ${fc}
                        fi
                    else
                        echo "${LOGPREFIX} Could not find ${dir}/${fc}/SampleSheet.csv."
                        launch_run ${dir} ${fc}
                    fi
                fi
            else
                echo "${LOGPREFIX} PASS: Sequencing not finished. Waiting for CopyComplete.txt"
            fi
        fi
        # else: ignore because format of folder name doesn't look like a run
    done
    echo "================================="
done
exit 0

# CopyComplete.txt ? --Y--> FastqComplete.txt ? --Y--> PASS
#       |                         |
#       N                         N
#       |                         |
#       V                         V
#      PASS                 Failed.txt ? ---------Y--> PASS
#                                 |
#                                 N
#                                 |
#                                 V
#                           LowPass.csv ? --------Y--> PASS
#                                 |
#                                 N
#                                 |
#                                 V
#                           SampleSheet.csv ? ----Y--> grep LowPass SampleSheet.csv ? ---------Y--> PASS
#                                 |                          |
#                                 |                    grep Cloud_Workflow SampleSheet.csv ? --Y--> PASS
#                                 |                          |
#                                 N                          N
#                                 |__________________________|
#                                               |
#                                               v
#                                           launch_run():
#                                           CopyComplete.txt ? --N--> PASS
#                                               |
#                                               Y
#                                               |
#                                               v
#                                 qsub dragen_bcl-convert_launcher.sh ${fc}
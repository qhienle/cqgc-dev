#!/bin/bash

# Watcher to detect and run pipeline for submitting cases to Emedgene.
# USAGE: bash emg_pipeline_watcher.sh
#        bash /mnt/spxp-app02/staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_pipeline_watcher.sh
#
# Logs for the process are written to these files:
#   - ${WORKDIR}/emg_watcher.log            : logs for this watcher
#   - ${WORKDIR}/${fc}/run_pipeline_emg.log : logs for each pipeline
#
# In addition, the following files mark the different steps of the pipeline:
#   - ${BASEDIR}/${FC}/CopyComplete.txt  : end of sequencing run
#   - ${BASEDIR}/${FC}/DemuxFailed.txt or "failed.txt" : sequencing has failed
#   - ${BASEDIR}/${FC}/DemuxStarted.txt  : bcl-convert is in progress
#   - ${WORKDIR}/${FC}/SampleSheet.csv   : bcl-convert is in progress
#   - ${BASEDIR}/${FC}/FastqComplete.txt : of BCL-conversion

HISEQ_R='/mnt/spxp-app02/staging/hiseq_raw'
SOFTDIR="/mnt/spxp-app02/staging2/soft/CQGC-utils"
BASEDIR='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs'
WORKDIR='/mnt/vs_nas_chusj/CQGC_PROD/fastqs'
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/LH00207R" "${HISEQ_R}/LH00336" "${HISEQ_R}/A00977" "${HISEQ_R}/LH00207R")
LOGPREFIX="[emg-watcher]"
LOGFILE="${WORKDIR}/emg_watcher.log"
NAPTIME=900

umask 002
printf "\n\n######\n%s %s %s\n######\n\n" ${LOGPREFIX} $( date "+%F@%T" ) $0 | tee -a ${LOGFILE}

# Run the pipeline
#
run_pipeline_emg() {
    local fc=$1
    local project=$2
    log="${WORKDIR}/${fc}/run_pipeline_emg.log"
    if [[ ! -d "${WORKDIR}/${fc}" ]]; then mkdir ${WORKDIR}/${fc}; fi
    cd ${WORKDIR}/${fc}
    echo "${LOGPREFIX} Get list of samples for run ${fc}" | tee -a ${log}
    python ${SOFTDIR}/Analysis.pipeline_illumina/list_run_samples.py ${fc}
    echo "${LOGPREFIX} Waiting for sequencing and demux to finish" | tee -a ${log}
    until [ -f "${BASEDIR}/${FC}/FastqComplete.txt" ]; do
        sleep ${NAPTIME}
    done
    echo "${LOGPREFIX} Demux has completed" | tee -a ${log}
    echo "${LOGPREFIX} Uploading samples to BaseSpace" | tee -a ${log}
    ## TODO: move FASTQ files from 1.fastq/PROJECT_NAME to 1.fastq/ before uploading when instrument is NovaSeq6000
    python ${SOFTDIR}/Analysis.pipeline_illumina/emg_upload_fastqs.py
    touch ${WORKDIR}/${FC}/UploadBsComplete.txt
    sleep ${NAPTIME}
    python ${SOFTDIR}/Analysis.pipeline_illumina/emg_make_batch.py --project ${project} >> ${WORKDIR}/${FC}/emg_make_batch.log 2>&1
}

for dir in ${WATCHDIRS[@]}; do
    for fc in $( ls ${dir} ); do
        # "RunName" in NovaSeq600's SampleSheet doesn't contain experiment name
        # (PRAG, Q1K,...) so better check Project Name in BaseSpace (with `bs`)
        #
        xp=$( bs -c cac1 list runs --format csv | grep ${fc} | cut -d, -f3 )
        if [[ ! -z ${xp} ]]; then
            if [[ -f "${WORKDIR}/${fc}/run_pipeline_emg.log" ]]; then
                echo "${LOGPREFIX} PASS: Pipeline already processed or is running for ${fc} | ${xp}" | tee -a ${LOGFILE}
            elif [[ ! -z $( echo ${xp} | grep 'PRAG' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for PRAG ${fc} | ${xp}" | tee -a ${LOGFILE}
                (run_pipeline_emg ${fc} 'prag') &
            elif [[ ! -z $( echo ${xp} | grep 'Q1K' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for Q1K ${fc} | ${xp}" | tee -a ${LOGFILE}
                (run_pipeline_emg ${fc} 'q1k') &
            elif [[ ! -z $( echo ${xp} | grep 'AOH' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for AOH ${fc} | ${xp}" | tee -a ${LOGFILE}
                (run_pipeline_emg ${fc} 'aoh') &
            elif [[ ! -z $( echo ${xp} | grep 'C4R' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for C4R ${fc} | ${xp}" | tee -a ${LOGFILE}
                (run_pipeline_emg ${fc} 'c4r') &
            else
                echo "${LOGPREFIX} Nothing to do for ${fc} (${xp})" | tee -a ${LOGFILE}
            fi
        else
            echo "${LOGPREFIX} Run ${fc} not found on BaseSpace." | tee -a ${LOGFILE}
        fi
    done
done
exit 0
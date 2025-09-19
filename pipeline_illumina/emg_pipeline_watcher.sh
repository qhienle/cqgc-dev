#!/bin/bash

# Watcher to detect and run pipeline for submitting cases to Emedgene.
# USAGE: bash run_pipeline_emg.sh
#        bash run_pipeline_emg.sh | tee -a /mnt/vs_nas_chusj/CQGC_PROD/fastqs/emg_watcher.log
#
# Logs for the process are written to these files:
#   - ${WORKDIR}/${fc}/run_pipeline_emg.log : logs for the pipeline
#   - ${WORKDIR}/emg_watcher.log            : logs for this watcher
#
# The following files mark the different steps of the pipeline:
#   - ${BASEDIR}/${FC}/CopyComplete.txt  : end of sequencing run
#   - ${BASEDIR}/${FC}/DemuxFailed.txt or "failed.txt" : sequencing has failed
#   - ${BASEDIR}/${FC}/DemuxStarted.txt  : bcl-convert is in progress
#   - ${WORKDIR}/${FC}/SampleSheet.csv   : bcl-convert is in progress
#   - ${BASEDIR}/${FC}/FastqComplete.txt : of BCL-conversion

HISEQ_R='/mnt/spxp-app02/staging/hiseq_raw'
BASEDIR='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs'
WORKDIR='/mnt/vs_nas_chusj/CQGC_PROD/fastqs'
SOFTDIR="/staging2/soft/CQGC-utils"
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/LH00207R" "${HISEQ_R}/LH00336" "${HISEQ_R}/A00977" "${HISEQ_R}/LH00207R")
LOGPREFIX="[emg-watcher]"
LOGFILE="${WORKDIR}/emg_watcher.log"
NAPTIME=900

umask 002
printf "\n\n######\n%s %s %s\n######\n\n" ${LOGPREFIX} $( date "+%F@%T" ) $0 #| tee -a ${LOGFILE}/

# Run the pipeline
#
run_pipeline_emg() {
    local fc=$1
    local project=$2
    log="${WORKDIR}/${fc}/run_pipeline_emg.log"
    if [[ ! -d "${WORKDIR}/${fc}" ]]; then mkdir ${WORKDIR}/${fc}; fi
    cd ${WORKDIR}/${fc}
    echo "Get list of samples for run ${fc}"
    python ${SOFTDIR}/Analysis.pipeline_illumina/list_run_samples.py ${fc}
    echo "${LOGPREFIX} Waiting for sequencing and demux to finish"
    # until [ -f "${BASEDIR}/${FC}/FastqComplete.txt" ]; do
    #     sleep ${NAPTIME}
    # done
    # echo "${LOGPREFIX} Demux has completed"
    # echo "Uploading samples to BaseSpace"
    # ## TODO: move FASTQ files from 1.fastq/PROJECT_NAME to 1.fastq/ before uploading when instrument is NovaSeq6000
    # python ${SOFTDIR}/Analysis.pipeline_illumina/emg_upload_fastqs.py
    # touch ${WORKDIR}/${FC}/UploadBsComplete.txt
    # sleep ${NAPTIME}
    # python ${SOFTDIR}/Analysis.pipeline_illumina/emg_make_batch.py >> ${WORKDIR}/${FC}/emg_make_batch.log 2>&1
}

for dir in ${WATCHDIRS[@]}; do
    for fc in $( ls ${dir} ); do
        #
        # "RunName" in NovaSeq600's SampleSheet don't contain experiment name
        # (PRAG, Q1K,...) so better check Project Name in BaseSpace (with `bs`)
        #
        xp=$( bs -c cac1 list runs --format csv | grep ${fc} | cut -d, -f3)
        if [[ ! -z ${xp} ]]; then
            if [[ -f "${WORKDIR}/${fc}/run_pipeline_emg.log" ]]; then
                echo "${LOGPREFIX} PASS: Pipeline already processed or is running for ${fc} | ${xp}"
            elif [[ ! -z $( echo ${xp} | grep 'PRAG' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for PRAG ${fc} | ${xp}"
                run_pipeline_emg ${fc} 'prag'
            elif [[ ! -z $( echo ${xp} | grep 'Q1K' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for Q1K ${fc} | ${xp}"
                run_pipeline_emg ${fc} 'q1k'
            elif [[ ! -z $( echo ${xp} | grep 'AOH' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for AOH ${fc} | ${xp}"
                run_pipeline_emg ${fc} 'aoh'
            elif [[ ! -z $( echo ${xp} | grep 'C4R' ) ]]; then
                echo "${LOGPREFIX} RUN: pipeline for C4R ${fc} | ${xp}"
                run_pipeline_emg ${fc} 'c4r'
            else
                true # echo "${LOGPREFIX} Nothing to do for ${fc} (${xp})"
            fi
        else
            true #echo "${LOGPREFIX} Run ${fc} not found on BaseSpace."
        fi
    done
done
exit 0
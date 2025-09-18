#!/bin/bash

# Watcher to detect and run pipeline for submitting cases to Emedgene.
# USAGE: bash run_pipeline_emg.sh
#        bash run_pipeline_emg.sh | tee -a /mnt/vs_nas_chusj/CQGC_PROD/fastqs/emg_watcher.log
#

HISEQ_R='/mnt/spxp-app02/staging/hiseq_raw'
BASEDIR='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs'
WORKDIR='/mnt/vs_nas_chusj/CQGC_PROD/fastqs'
SOFTDIR="/staging2/soft/CQGC-utils"
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/LH00207R" "${HISEQ_R}/LH00336" "${HISEQ_R}/A00977" "${HISEQ_R}/LH00207R")

LOGPREFIX="[emg-watcher]"
LOGFILE="${WORKDIR}/emg_watcher.log"

# Run the pipeline
#
run_pipeline_emg() {
    local fc=$1
    local project=$2
    cd ${WORKDIR}/${fc}
    echo "${LOGPREFIX} ${project} bash ${SOFTDIR}/Analysis.pipeline_illumina/run_pipeline_prag.sh ${fc} 2>&1 | tee ${WORKDIR}/${fc}/run_pipeline_prag.log"
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
                echo "${LOGPREFIX} PASS: Pipeline already processed or is running"
            elif [[ ! -z $( echo ${xp} | grep 'PRAG' ) ]]; then
                echo "${LOGPREFIX} Launching pipeline for PRAG ${fc} | (${xp})"
                run_pipeline_emg ${fc} 'prag'
            elif [[ ! -z $( echo ${xp} | grep 'Q1K' ) ]]; then
                echo "${LOGPREFIX} Launching pipeline for Q1K ${fc} | (${xp})"
                run_pipeline_emg ${fc} 'q1k'
            elif [[ ! -z $( echo ${xp} | grep 'AOH' ) ]]; then
                echo "${LOGPREFIX} Launching pipeline for AOH ${fc} | (${xp})"
                run_pipeline_emg ${fc} 'aoh'
            elif [[ ! -z $( echo ${xp} | grep 'C4R' ) ]]; then
                echo "${LOGPREFIX} Launching pipeline for C4R ${fc} | (${xp})"
                run_pipeline_emg ${fc} 'c4r'
            else
                echo "${LOGPREFIX} Nothing to do for ${fc} (${xp})"
            fi
        else
            echo "${LOGPREFIX} Run ${fc} not found on BaseSpace."
        fi
    done
done

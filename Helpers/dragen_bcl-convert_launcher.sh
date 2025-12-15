#!/bin/bash
#$ -N 'dragen_bcl_convert'
#$ -j y
#$ -l dragen=1

# Submit BCL-convert job to DRAGEN servers, from spxp-app02.
# qsub /staging2/soft/CQGC-utils/Helpers/dragen_bcl_convert_launcher.sh <FLOWCELL>
# qsub \
#   -v BASEDIR="/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/LH00336
#   -v WORKDIR="/mnt/vs_nas_chusj/CQGC_PROD/fastqs" \
#   /staging2/soft/CQGC-utils/Helpers/dragen_bcl_convert_launcher.sh <FLOWCELL>
# or export environment variables:
#   export BASEDIR="/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/LH00336
#   export WORKDIR="/mnt/vs_nas_chusj/CQGC_PROD/fastqs"

# <FLOWCELL> argument is REQUIRED
# DemuxStarted.txt file in ${BASEDIR}/${FC}/ marks a bcl-convert in progress
# FastqComplete.txt file created by DRAGEN marks end of process
# DemuxFailed.txt file in ${BASEDIR}/${FC}/ marks a bcl-convert failure

umask 002

mark_failed() {
    # In case of error, mark run with "DemuxFailed.txt"
    echo "$1" >&2
    touch DemuxFailed.txt
    exit 1
}

if [[ -z ${1} ]]; then
    mark_failed "ERROR: Flowcell or run name not provided!"
else
    FC=${1}
    a=($(echo ${FC} | tr '_' '\n'))
    FC_SHORT="${a[1]}_${a[2]}"
fi

# Set BASEDIR and WORKDIR if environment variables not exported
if [[ -z ${BASEDIR} ]]; then 
    if [[ -f "/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/${a[1]}/${FC}/CopyComplete.txt" ]]; then
        BASEDIR="/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/${a[1]}"
    elif [[ -f "/mnt/spxp-app02/staging/hiseq_raw/${a[1]}/${FC}/CopyComplete.txt" ]]; then
        BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/${a[1]}"
    else
        mark_failed "ERROR: Could not find: '${BASEDIR}/${FC}/CopyComplete.txt'"
    fi
else
    if [[ ! -f "${BASEDIR}/${FC}/CopyComplete.txt" ]]; then
        mark_failed "ERROR: Could not find '${FC}/CopyComplete.txt' in specified BASEDIR '${BASEDIR}'"
    fi
fi
if [[ -z ${WORKDIR} ]]; then
    # WORKDIR="/mnt/spxp-app02/staging2/dragen"
    WORKDIR="/mnt/vs_nas_chusj/CQGC_PROD/fastqs"
fi
# Set OUTDIR and check that there isn't a prior demux in progress
OUTDIR="${WORKDIR}/${FC}/1.fastq"


# Set SampleSheet.csv. Fix funky characters in SampleNames.txt and other Nanuq files
# TODO: Fix this in get_nanuq_files.py
# Run bcl-convert depending on the instrument and SampleSheet
samplesheet="${WORKDIR}/${FC}/SampleSheet.csv"
if [[ -f ${samplesheet} ]]; then
    touch ${BASEDIR}/${FC}/DemuxStarted.txt
    if [[ "${FC_SHORT}" =~ ^A00* ]]; then
        echo "Run dragen BCL-convert for ${FC}"
        dragen \
            --bcl-conversion-only true \
            --bcl-input-directory ${BASEDIR}/${FC} \
            --output-directory ${OUTDIR} \
            --sample-sheet ${samplesheet} \
            --bcl-only-matched-reads true \
            --bcl-sampleproject-subdirectories true \
            >> ${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1
    elif [[ "${FC_SHORT}" =~ ^LH00* ]]; then
        echo "Run dragen BCL-convert for ${FC}"
        dragen \
            --bcl-conversion-only true \
            --bcl-input-directory ${BASEDIR}/${FC} \
            --output-directory ${OUTDIR} \
            --sample-sheet ${samplesheet} \
            --bcl-only-matched-reads true \
            >> ${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1
    else
        mark_failed "ERROR: Could not determine instrument series for ${FC}"
    fi
    if [[ $? -eq 0 ]]; then
        cp ${OUTDIR}/Logs/FastqComplete.txt ${BASEDIR}/${FC}
        mv ${OUTDIR}/streaming_log_${USER}.csv ${OUTDIR}/Logs/
        mv ${OUTDIR}/dragen* ${OUTDIR}/Logs/
    else
        mark_failed "ERROR: Something went wrong! DRAGEN exit code not equal to zero."
        touch ${BASEDIR}/${FC}/FastqComplete.txt # or FastqFailed.txt ?
    fi
else
    mark_failed "ERROR: Did you forget to setup environment variables or provide the SampleSheet?\n" >&2
fi

exit 0

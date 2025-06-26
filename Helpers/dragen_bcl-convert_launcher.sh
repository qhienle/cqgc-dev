#!/bin/bash
#$ -N 'dragen_bcl_convert'
#$ -j y
#$ -l dragen=1

# Submit BCL-convert job to DRAGEN servers, from spxp-app02.
# qsub /staging2/soft/CQGC-utils/Helpers/dragen_bcl_convert.sh <FLOWCELL>
# qsub \
#   -v BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/LH00336
#   -v WORKDIR="/mnt/spxp-app02/staging2/dragen" \
#   /staging2/soft/CQGC-utils/Helpers/dragen_bcl_convert.sh <FLOWCELL>
# or export environment variables:
#   export BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/LH00336
#   export WORKDIR="/mnt/spxp-app02/staging2/dragen"

# <FLOWCELL> argument is REQUIRED
# DemuxStarted.txt file in ${BASEDIR}/${FC}/ marks a bcl-convert in progress
# FastqComplete.txt file created by DRAGEN marks end of process

umask 002

if [[ -z ${1} ]]; then
    echo "ERROR: Flowcell or run name not provided!" >&2
    exit 1
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
        echo "ERROR: Could not find: '${BASEDIR}/${FC}/CopyComplete.txt'"
        exit 1
    fi
else
    if [[ ! -f "${BASEDIR}/${FC}/CopyComplete.txt" ]]; then
        echo "ERROR: Could not find '${FC}/CopyComplete.txt' in specified BASEDIR '${BASEDIR}'"
        exit 1
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
        echo "ERROR: Could not determine instrument series for ${FC}" >&2
        exit 1
    fi
    echo "TEST exit code is $?"
    cp ${OUTDIR}/Logs/FastqComplete.txt ${BASEDIR}/${FC}
    mv ${OUTDIR}/streaming_log_${USER}.csv ${OUTDIR}/Logs/
    mv ${OUTDIR}/dragen* ${OUTDIR}/Logs/
else
    echo "ERROR: Did you forget to setup environment variables or provide the SampleSheet?\n" >&2
    exit 1
fi

exit 0

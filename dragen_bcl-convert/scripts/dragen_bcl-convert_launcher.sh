#!/bin/bash
#$ -N 'dragen_bcl_convert'
#$ -j y
#$ -l dragen=1

# Submit BCL-convert by DRAGEN on spxp-app04, from spxp-app02.
# qsub /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl_convert.sh <FLOWCELL>

if [ -z ${1} ]; then
    echo "ERROR: Flowcell or run name not provided!" >&2
    exit 1
fi

FC=${1}
a=($(echo ${FC} | tr '_' '\n'))
FC_SHORT="${a[1]}_${a[2]}"
BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/${a[1]}"
WORKDIR="/mnt/spxp-app02/staging2/dragen"
OUTDIR="${WORKDIR}/${FC}/1.fastq"

if [ ! -d ${WORKDIR}/${FC} ]; then 
    mkdir ${WORKDIR}/${FC}
fi

if [[ "${FC_SHORT}" =~ ^A00* ]]; then
    echo "Get Nanuq files for NoveSeq6000 run ${FC}"
    echo "Run dragen BCL-convert for ${FC}"
    dragen \
        --bcl-conversion-only true \
        --bcl-input-directory ${BASEDIR}/${FC} \
        --output-directory ${OUTDIR} \
        --sample-sheet ${WORKDIR}/${FC}/SampleSheet.csv \
        --bcl-only-matched-reads true \
        --bcl-sampleproject-subdirectories true \
        --force \
        >> ${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1
elif [[ "${FC_SHORT}" =~ ^LH00* ]]; then
    echo "Get Nanuq files for NoveSeqXPlus run ${FC}"
    echo "Run dragen BCL-convert for ${FC}"
    dragen \
        --bcl-conversion-only true \
        --bcl-input-directory ${BASEDIR}/${FC} \
        --output-directory ${OUTDIR} \
        --sample-sheet ${WORKDIR}/${FC}/SampleSheet.csv \
        --bcl-only-matched-reads true \
        --force \
        >> ${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1
else
    echo "ERROR: Could not determine instrument series for ${FC}" >&2
    exit 1
fi
touch ${WORKDIR}/${FC}/DemuxComplete.txt
mv ${OUTDIR}/streaming_log_${USER}.csv ${OUTDIR}/Logs/
mv ${OUTDIR}/dragen* ${OUTDIR}/Logs/

exit 0

#!/bin/bash
#$ -N 'dragen_bcl_convert'
#$ -j y
#$ -l dragen=1

# Submit BCL-convert by DRAGEN on spxp-app04, from spxp-app02.
# qsub /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl_convert.sh <FLOWCELL>

FC=$1
if [ -z ${FC} ]; then
    echo "ERROR: Flowcell or run name not provided!"
    exit
else
    a=($(echo ${FC} | tr '_' '\n'))
    FC_SHORT="${a[1]}_${a[2]}"
    BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/${a[1]}"
    WORKDIR="/mnt/spxp-app02/staging2/dragen"
    OUTDIR="${WORKDIR}/${FC}/1.fastq"
    if [[ "${FC_SHORT}" =~ ^A00* ]]; then
        INSTR="NoveSeq6000"
        echo "Get Nanuq files for ${INSTR} run ${FC}"
        python /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py --run ${FC_SHORT}
        echo "Run dragen BCL-convert for ${FC}"
        # dragen \
        #     --bcl-conversion-only true \
        #     --bcl-input-directory ${BASEDIR}/${FC} \
        #     --output-directory ${OUTDIR} \
        #     --sample-sheet ${WORKDIR}/${FC}/SampleSheet.csv \
        #     --bcl-only-matched-reads true \
        #     --bcl-sampleproject-subdirectories true \
        #     --force \
        #     >> ${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1
    elif [[ "${FC_SHORT}" =~ ^LH00* ]]; then
        INSTR="NoveSeqXPlus"
        echo "Get Nanuq files for ${INSTR} run ${FC}"
        python /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py --orient-index2 --run ${FC_SHORT}
        echo "Run dragen BCL-convert for ${FC}"
        # dragen \
        #     --bcl-conversion-only true \
        #     --bcl-input-directory ${BASEDIR}/${FC} \
        #     --output-directory ${OUTDIR} \
        #     --sample-sheet ${WORKDIR}/${FC}/SampleSheet.csv \
        #     --bcl-only-matched-reads true \
        #     --force \
        #     >> ${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1
    else
        echo "ERROR: Could not determine instrument series for ${FC}"
        exit
    fi
    mv ${OUTDIR}/streaming_log_${USER}.csv ${OUTDIR}/Logs/
    mv ${OUTDIR}/dragen* ${OUTDIR}/Logs/
fi
exit

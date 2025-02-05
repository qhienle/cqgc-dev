#!/bin/bash

# Run pipeline for PRAGMatIQ
# USAGE: bash run_pipeline_prag.sh <RUN>
#        bash run_pipeline_prag.sh 20241220_LH00336_0145_B22MKV5LT3
#        bash run_pipeline_prag.sh ${FC} >>${WORKDIR}/${FC}/${FC_SHORT}.bcl-convert.log 2>&1

## 0. Mise en place de l'environnement de travail
FC=${1}
a=($(echo ${FC} | tr '_' '\n'))
FC_SHORT="${a[1]}_${a[2]}"
BASEDIR="/mnt/spxp-app02/staging/hiseq_raw/${a[1]}"
WORKDIR="/mnt/spxp-app02/staging2/dragen"
SOFTDIR="/staging2/soft/CQGC-utils"
NAPTIME=900

if [ ! -d ${WORKDIR}/${FC} ]; then
    mkdir ${WORKDIR}/${FC}
fi
cd ${WORKDIR}/${FC}

## 1. Collecter les informations sur les familles
echo "Get list of samples for run ${FC}"
python ${SOFTDIR}/Analysis.pipeline_illumina/list_run_samples.py ${FC}


## 2. Préparer les FASTQs
## 2.1. La déconvolution et conversion des BCLs en FASTQs devrait se faire automatiquement
# echo "Waiting for ${BASEDIR}/${FC}/CopyComplete.txt"
# until [ -f ${BASEDIR}/${FC}/CopyComplete.txt ]
# do
#     printf '.'
#     sleep ${NAPTIME}
# done
# echo
# echo "Sequencing has completed."
# echo "Getting SampleSheet from Nanuq"
# python ${SOFTDIR}/Helpers/get_nanuq_files.py --run ${FC_SHORT} 
# qsub -V ${SOFTDIR}//Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_launcher.sh ${FC}

### 2.2. Téléverser les FASTQs sur BaseSpace
until [ -f ${WORKDIR}/${FC}/DemuxComplete.txt ]
do
    printf '.'
    sleep ${NAPTIME}
done
echo
echo "Demux has completed. Uploading samples to BaseSpace"
## TODO: move FASTQ files from 1.fastq/PROJECT_NAME to 1.fastq/ before uploading when instrument is NovaSeq6000
python ${SOFTDIR}/Analysis.pipeline_illumina/emg_upload_fastqs.py
touch ${WORKDIR}/${FC}/UploadBsComplete.txt


## 3. Créer les cas sur Emedgene 
###  3.1. Générer le fichier "emg_batch_manifest.csv" `emg_make_batch_from_nanuq.py ${FC}`
###  3.2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
# TODO: use `bs` to check if files are on BaseSpace for each sample
# python ${SOFTDIR}/Analysis.pipeline_illumina/emg_make_batch_from_nanuq.py ${FC_SHORT} >> ${WORKDIR}/${FC}/emg_make_batch.log 2>&1


## 4. Collecter les metriques
# python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_collect_samples_metrics.py ${FC}


## 5. Archiver les résultats
# ssh hien@narval.calculquebec.ca

# samples=()
# total_samples=${#samples[@]}

# cd ${HOME}/projects/ctb-rallard/COMMUN/PRAGMatIQ-EMG
# cp archive_PRAGMatIQ.log archive_PRAGMatIQ.log0

# count=1
# for sample in ${samples[@]}; do
#     echo "Processing sample ${count}/${total_samples}, ${sample}"
#     aws s3 --profile emedgene cp s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/${sample}/ ./archives/${sample} --recursive
#     echo ${sample} $( date +'%Y-%m-%d %T' ) >> archive_PRAGMatIQ.log
#     ((count++))
# done
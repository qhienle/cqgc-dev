#!/bin/bash

# Run pipeline for PRAGMatIQ
# USAGE: bash run_pipeline_prag.sh <RUN>
#        bash run_pipeline_prag.sh 20241220_LH00336_0145_B22MKV5LT3
#        bash run_pipeline_prag.sh ${FC} >>${WORKDIR}/${FC}/${FC}.pipeline_prag.log 2>&1

## 0. Mise en place de l'environnement de travail
SOFTDIR="/staging2/soft/CQGC-utils"
NAPTIME=900
FC=${1}
a=($(echo ${FC} | tr '_' '\n'))
if [[ -z ${BASEDIR} ]]; then
    BASEDIR="/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/${a[1]}"
else
    BASEDIR=${BASEDIR} # export BASEDIR="/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/${a[1]}"
fi
if [[ -z ${WORKDIR} ]]; then
    WORKDIR="/mnt/vs_nas_chusj/CQGC_PROD/fastqs"
else
    WORKDIR=${WORKDIR} # export WORKDIR="/mnt/vs_nas_chusj/CQGC_PROD/fastqs"
fi

## 1. Collecter les informations sur les familles dans samples_list.csv,
## fournit aussi la liste des échantillons à téléverser sur BaseSpace (étape 3)
echo "Get list of samples for run ${FC}"
if [[ ! -d "${WORKDIR}/${FC}" ]]; then
    mkdir ${WORKDIR}/${FC};
fi
cd ${WORKDIR}/${FC}
python ${SOFTDIR}/Analysis.pipeline_illumina/list_run_samples.py ${FC}


## 2. Déconvolution et conversion des BCLs en FASTQs: devrait se faire 
## automatiquement par dragen_bcl-convert_watcher.sh
echo "Waiting for sequencing to finish"
until [ -f ${BASEDIR}/${FC}/CopyComplete.txt ]
do
    printf '.'
    sleep ${NAPTIME}
done
echo "Sequencing has completed. Waiting for Demux to finish"
until [ -f "${BASEDIR}/${FC}/FastqComplete.txt" ]
do
    printf '.'
    sleep ${NAPTIME}
done
echo "Demux has completed"


## 3. Téléverser les FASTQs sur BaseSpace
echo "Uploading samples to BaseSpace"
python ${SOFTDIR}/Analysis.pipeline_illumina/emg_upload_fastqs.py
touch ${WORKDIR}/${FC}/UploadBsComplete.txt


## 4. Créer les cas sur Emedgene 
###  4.1. Générer le fichier "emg_batch_manifest.csv" `emg_make_batch_from_nanuq.py ${FC}`
###  4.2. Glisser-déposer dans Emedgene le fichier "emg_batch_manifest.csv"
# TODO: use `bs` to check if files are on BaseSpace for each sample
# Fails if launched immediately after UploadBsComplete.
# Wait a while for the last upload to register with BaseSpace.
sleep ${NAPTIME}
python ${SOFTDIR}/Analysis.pipeline_illumina/emg_make_batch.py >> ${WORKDIR}/${FC}/emg_make_batch.log 2>&1


## 5. Collecter les metriques
# python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/emg_collect_samples_metrics.py ${FC}


## 6. Archiver les résultats
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
#!/usr/bin/bash
#
# USAGE: archive_PRAGMatIQ.sh sample1 sample2 sample3 ...
#
# Archive PRAGMatIQ Runs from Emedgenes AWS S3 bucket to narval
# ssh narval.calculquebec.ca
#
samples=$@
total_samples=${#@}
samples=(MO-24-000214 MO-24-000215 24-00002-T1 24-00087-T1 24-00089-T1)
total_samples=${#samples[@]}

cd ${HOME}/projects/ctb-rallard/COMMUN/PRAGMatIQ-EMG
cp archive_PRAGMatIQ.log archive_PRAGMatIQ.log0

count=1
for sample in ${samples[@]}; do
    echo "Processing sample ${count}/${total_samples}, ${sample}"
    aws s3 --profile emedgene cp s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/${sample}/ ./archives/${sample} --recursive
    echo ${sample} $( date +'%Y-%m-%d %T' ) >> archive_PRAGMatIQ.log
    ((count++))
done

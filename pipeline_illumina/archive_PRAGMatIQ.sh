#!/usr/bin/bash
#
# USAGE: archive_PRAGMatIQ.sh sample1 sample2 sample3 ...
#
# Archive PRAGMatIQ Runs from Emedgenes AWS S3 bucket to narval
# ssh narval.calculquebec.ca
#
cd ${HOME}/projects/ctb-rallard/COMMUN/PRAGMatIQ-EMG
count=1
total_samples=${#@}
for sample in ${@}; do
    echo "Processing sample ${count}/${total_samples}, ${sample}"
    aws s3 --profile emedgene cp s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/${sample}/ ./archives/${sample} --recursive
    echo ${sample} $( date +'%Y-%m-%d %T' ) >> archive_PRAGMatIQ.log
    ((count++))
done

#!/usr/bin/bash
#
# USAGE: archive_PRAGMatIQ.sh sample1 sample2 sample3 ...
#
# Archive PRAGMatIQ Runs from Emedgenes AWS S3 bucket to narval
# ssh narval.calculquebec.ca
#
samples=$@
total_samples=${#@}
samples=(GM232116 GM232117 GM232375 GM232377 GM232376 3103495962 3103495961 3113525391 3113525505 3113525463 3113545412 3113545511 23-09628-T1 23-09564-T1 23-09670-T1)
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

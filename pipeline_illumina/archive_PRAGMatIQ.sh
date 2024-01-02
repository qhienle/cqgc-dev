#!/usr/bin/bash
#
# USAGE: archive_PRAGMatIQ.sh sample1 sample2 sample3 ...
#
# Archive PRAGMatIQ Runs from Emedgenes AWS S3 bucket to narval
# ssh narval.calculquebec.ca
#
samples=$@
total_samples=${#@}
#samples=(GM232704 GM232705 GM232708 GM232658 GM232650 GM232675 GM232700 GM232674 GM232673 GM232723 GM232733 GM232732 3123698762 3123698754 3123698744 23-10447-T1 23-10583-T1 23-10591-T1)
samples=(GM232650 GM232675 GM232700 GM232674 GM232673 GM232723 GM232733 GM232732 3123698762 3123698754 3123698744 23-10447-T1 23-10583-T1 23-10591-T1)
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

#!/usr/bin/bash
#
# USAGE: archive_PRAGMatIQ.sh sample1 sample2 sample3 ...
#
# Archive PRAGMatIQ Runs from Emedgenes AWS S3 bucket to narval
# ssh narval.calculquebec.ca
#
samples=$@
total_samples=${#@}
samples=(GM240041 GM240042 GM240040 GM240068 GM240069 GM240070 GM240123 GM240125 GM240126 4013833470 4013850929 4013849723 MO-23-014288 MO-24-000888 MO-24-000887 24-00666-T1 24-00667-T1 24-00663-T1)
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

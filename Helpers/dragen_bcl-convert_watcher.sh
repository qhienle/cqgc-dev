#!/bin/bash

# Watch for new sequencing runs to launch DRAGEN BCL-Convert
# USAGE: Schedule launch with `crontab -e`, as `cqgc_user` (see below)
#        bash dragen_bcl-convert_watcher.sh | tee -a ${LOGFILE}
#        bash /staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_watcher.sh | tee -a /mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/dragen_bcl-convert_watcher.log

# # As "cqgc_user", run dragen_bcl-convert_watcher.sh and rotate its' log every Sunday
# # sudo su cqgc_user --login && crontab -l|-e
# */30 * * * * /usr/bin/bash -l -c "/staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_watcher.sh >> /mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/dragen_bcl-convert_watcher.log 2>&1"
# 59 23 * * Sun path='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs/dragen_bcl-convert_watcher'; flog="${path}.log"; dlog="${path}-$( date +"\%Y\%m\%d" ).log"; mv ${flog} ${dlog}; gunzip ${flog}.tar.gz; tar -rf ${flog}.tar ${dlog}; gzip ${flog}.tar; rm -f ${dlog};

# Scan BCL output dirs (BASEDIR) for new runs (FC) to demux
# Check SampleSheet, launch demux if not LowPass or Cloud_Workflow
#   Get SampleSheet if it doesn't exist and launch demux
# Demux outputs are written to ${WORKDIR}/dragen/${FC}

# File ${BASEDIR}/${FC}/CopyComplete.txt marks end of sequencing run
# File ${BASEDIR}/${FC}/DemuxFailed.txt or "failed.txt" marks that sequencing has failed
# File ${BASEDIR}/${FC}/DemuxStarted.txt marks that bcl-convert is in progress
#   DemuxStarted.txt is written to ${BASEDIR} by dragen_bcl-convert_launcher.sh
# File ${WORKDIR}/${FC}/SampleSheet.csv also marks that bcl-convert is in progress
# File ${BASEDIR}/${FC}/FastqComplete.txt marks end of BCL-conversion
#   FastqComplete.txt is copied to ${BASEDIR} by dragen_bcl-convert_launcher.sh

# DEPENDENCIES:
#   /staging2/soft/CQGC-utils/Helpers/dragen_bcl-convert_launcher.sh
#   /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py

SOFTDIR="/mnt/spxp-app02/staging2/soft/CQGC-utils"
BASEDIR='/mnt/vs_nas_chusj/CQGC_PROD/sequenceurs'
HISEQ_R='/mnt/spxp-app02/staging/hiseq_raw'
WORKDIR='/mnt/vs_nas_chusj/CQGC_PROD/fastqs'
LOGFILE="${WORKDIR}/dragen_bcl-convert_watcher.log"
LOGPREFIX="[bcl-watcher]"
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/LH00207R" "${HISEQ_R}/LH00336" "${HISEQ_R}/A00977" "${HISEQ_R}/LH00207R")

umask 002
printf "\n\n######\n%s %s %s\n######\n\n" ${LOGPREFIX} $( date "+%F@%T" ) $0 #| tee -a ${LOGFILE}/

launch_run() {
    # Run dragen_bcl-convert_launcher.sh and gather Demultiplex_Stats when done
    # FastqComplete.txt copied to ${BASEDIR} by dragen_bcl-convert_launcher.sh
    # and marks the run as done.
    local dir="$1"
    local fc="$2"
    
    # Setting up environments for qsub and conda
    . /mnt/spxp-app02/staging2/soft/GE2011.11p1/SGE_ROOT/default/common/settings.sh
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate CQGC-utils
    cd ${WORKDIR}/${fc}
    echo "${LOGPREFIX} ${fc} RUN: Launching BCL-convert with qsub..."
    # echo "echo 'qsub moot launcher'" | qsub -V -o "${WORKDIR}/${fc}/qsub_out.txt" -e "${WORKDIR}/${fc}/qsub_err.txt" # for testing
    qsub -V -o "${WORKDIR}/${fc}/qsub_out.txt" -e "${WORKDIR}/${fc}/qsub_err.txt" ${SOFTDIR}/Helpers/dragen_bcl-convert_launcher.sh ${fc}
    touch ${dir}/${fc}/DemuxStarted.txt
    echo "Waiting for Demux to finish" >>qsub_out.txt 2>&1
    until [ -f "${dir}/${fc}/FastqComplete.txt" ]; do
        printf '.' >>qsub_out.txt 2>&1
        sleep 900 # 15 minutes
    done
    echo -e "\nDemux has completed. Gathering demultiplexing statstics for QC..." >>qsub_out.txt 2>&1
    bash ${SOFTDIR}/Analysis.dragen_bcl-convert/scripts/cp_RunInfo_Stats.sh ${dir} ${WORKDIR} ${fc} >>qsub_out.txt 2>&1
    python ${SOFTDIR}/Analysis.dragen_bcl-convert/scripts/qc_demultiplex_stats.py --file ${WORKDIR}/${fc}/1.fastq/Reports/Demultiplex_Stats.csv >>qsub_out.txt 2>&1
    ${SOFTDIR}/QC.dragen_demultiplexing/index_stats.py -d 1.fastq/Reports > 1.fastq/Reports/${fc}.index_stats.out # $QC_DIR="4.qualite/"
    conda deactivate
}

for dir in ${WATCHDIRS[@]}; do
    # echo "${LOGPREFIX} ====== Scanning ${dir}... ======"
    for fc in $( ls ${dir} ); do
        parts=($(echo ${fc} | tr '_' '\n'))
        if [[ ${#parts[@]} -eq 4 ]]; then
            # `bs` for debug. Probably doesn't work in crontab
            echo "${LOGPREFIX} $( bs -c cac1 list runs --format csv | grep ${fc} | sed -e 's/,/ \| /g')" 
            # echo "${LOGPREFIX} ${fc}"
            if [[ -f "${dir}/${fc}/CopyComplete.txt" ]]; then
                # echo "${LOGPREFIX} CopyComplete.txt indicates that sequencing has finished"
                # Check if bcl-convert needed (not previously demuxed, not failed, not LowPass)
                if [[ -f "${dir}/${fc}/FastqComplete.txt" ]]; then
                    echo "${LOGPREFIX} ${fc} PASS: FastqComplete.txt indicates that run has already been processed."
                elif [[ -f "${dir}/${fc}/DemuxFailed.txt" ]] ||  [[ -f "${dir}/${fc}/failed.txt" ]]; then
                    echo "${LOGPREFIX} ${fc} PASS: Failed.txt marks a failed Run."
                elif [[ -f "${dir}/${fc}/DemuxStarted.txt" ]]; then
                    echo "${LOGPREFIX} ${fc} PASS: DemuxStarted.txt marks a bcl-convert process in progress."
                elif [[ -f "${dir}/${fc}/LowPass*.csv" ]]; then
                    echo "${LOGPREFIX} ${fc} PASS: Found what looks like a LowPass SampleSheet."
                else
                    # Check SampleSheet, for LowPass or Cloud_Workflow.
                    if [[ -f "${dir}/${fc}/SampleSheet.csv" ]]; then
                        if grep -q "LowPass" "${dir}/${fc}/SampleSheet.csv"; then
                            echo "${LOGPREFIX} ${fc} PASS: SampleSheet for LowPass"
                        elif grep -q "Cloud_Workflow," "${dir}/${fc}/SampleSheet.csv"; then
                            echo "${LOGPREFIX} ${fc} PASS: SampleSheet indicates a Cloud_Workflow"
                            # TODO: Delete Run (if for TSO500)?
                        else
                            echo "${LOGPREFIX} ${fc} LAUNCH: SampleSheet found, not for LowPass or Cloud_Workflow."
                            launch_run ${dir} ${fc} &
                        fi
                    else
                        echo "${LOGPREFIX} ${fc} LAUNCH: Could not find ${dir}/${fc}/SampleSheet.csv. Getting files from Nanuq..."
                        mkdir ${WORKDIR}/${fc}
                        cd ${WORKDIR}/${fc}
                        python3 ${SOFTDIR}/Helpers/get_nanuq_files.py --run ${fc}
                        dos2unix ${WORKDIR}/${fc}/Sample*
                        launch_run ${dir} ${fc} &
                    fi
                fi
            else
                echo "${LOGPREFIX} ${fc} PASS: Sequencing not finished. Waiting for CopyComplete.txt"
            fi
        fi
        # else: ignore because format of folder name doesn't look like a run
    done
done
exit 0

# CopyComplete.txt ? --Y--> FastqComplete.txt ? --Y--> PASS
#       |                         |
#       N                         N
#       V                         V
#      PASS                 Failed.txt ? ---------Y--> PASS
#                                 |
#                                 N
#                                 V
#                         DemuxStarted.txt ? -----Y--> PASS
#                                 |
#                                 N
#                                 V
#                           LowPass.csv ? --------Y--> PASS
#                                 |
#                                 N
#                                 V
#                          SampleSheet.csv ? -----Y--> grep LowPass SampleSheet.csv ? ---------Y--> PASS
#                                 |                          |
#                                 N                          N
#                                 V                          V
#                          get_nanuq_files.py          grep Cloud_Workflow SampleSheet.csv ? --Y--> PASS
#                                 |                          |
#                                 |                          N
#                                 |__________________________|
#                                               V
#                                           launch_run(): qsub dragen_bcl-convert_launcher.sh ${fc}
#                                                                  V
#                                                         cp_RunInfo_Stats.sh
#                                                                  V
#                                                         qc_demultiplex_stats.py
#                                                                  V
#                                                         index_stats.py

# ```mermaid
# ---
# config:
#   theme: mc
#   look: classic
#   layout: dagre
# ---
# flowchart TD
#     A[Start: Scan WATCHDIRS for runs] --> B{Is folder name valid?}
#     B -->|No| C[fa:fa-ban PASS: Ignore folder]
#     B -->|Yes| D{CopyComplete.txt exists?}
#     D -->|No| E[fa:fa-ban PASS: Sequencing not finished]
#     D -->|Yes| F{FastqComplete.txt exists?}
#     F -->|Yes| G[fa:fa-ban PASS: Run already processed]
#     F -->|No| H{DemuxFailed.txt or failed.txt exists?}
#     H -->|Yes| I[fa:fa-ban PASS: Run marked as failed]
#     H -->|No| J{DemuxStarted.txt exists?}
#     J -->|Yes| K[fa:fa-ban PASS: Demux in progress]
#     J -->|No| L{LowPass.csv exists?}
#     L -->|Yes| M[fa:fa-ban PASS: LowPass SampleSheet found]
#     L -->|No| N{SampleSheet.csv exists?}
#     N -->|No| O[Fetch SampleSheet using get_nanuq_files.py]
#     O --> P[Convert SampleSheet to Unix format]
#     P --> Q[Launch dragen_bcl-convert_launcher.sh]
#     N -->|Yes| R{SampleSheet contains LowPass or Cloud_Workflow?}
#     R -->|Yes| S[fa:fa-ban PASS: LowPass or Cloud_Workflow SampleSheet]
#     R -->|No| Q
#     Q --> T[Wait for FastqComplete.txt]
#     T --> U[Gather demultiplexing stats]
#     U --> V[Run cp_RunInfo_Stats.sh]
#     V --> W[Run qc_demultiplex_stats.py]
#     W --> Y[Run index_stats.py]
#     C & E & G & I & K & M & S --> X[Continue scanning other folders]
#     Y --> X
# ```


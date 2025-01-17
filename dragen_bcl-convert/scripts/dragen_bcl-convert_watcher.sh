#!/bin/bash

# Watch for new sequencing runs to launch DRAGEN BCL-Convert and extract stats.
# This version to be launched using crontab, instead of systemd for the service
# version of this script `dragen_bcl-convert_watcherd.sh`.
# USAGE :
#       launch in crontab
#       nohup /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_watcherd.sh > /dev/null 2>&1 &

# Scan BCL output dirs (BASEDIR) for new runs (FC) to demux
# Skip runs for LowPass (check if SampleSheet exists (LowPass))
# File ${BASEDIR}/${FC}/CopyComplete.txt marks end of sequencing run
# Demux outputs are written to ${WORKDIR}/dragen/${FC}
# File ${WORKDIR}/dragen/${FC}/DemuxComplete.txt marks end of BCL-conversion


BASEDIR='/staging/hiseq_raw'
WORKDIR='/staging2/dragen'
WATCHDIRS=("${BASEDIR}/A00516" "${BASEDIR}/LH00336" "${BASEDIR}/A00977" "${BASEDIR}/LH00207R" "/mnt/vs_nas_chusj/SPXP_APP02_NFS/LH00336")
NAPTIME=900

launch_run() {
    # Check if sequencing is finished (CopyComplete.txt) and that run
    # is not already being processed by another instance of this script
    # which creates output dir and Sample* files through `get_nanuq_files.py``
    local dir="$1"
    local fc="$2"
    parts=($(echo ${fc} | tr '_' '\n'))
    fc_short="${parts[1]}_${parts[2]}"
    if [[ -f "${dir}/${fc}/CopyComplete.txt" ]]; then
        echo "Found file ${dir}/${fc}/CopyComplete.txt indicating that sequencing has finished"
        if [[ -d ${WORKDIR}/${fc} ]]; then
            echo "PASS: Run already processed in ${WORKDIR}/${fc}"
        else
            echo "Getting SampleSheet and other files from Nanuq..."
            python /staging2/soft/CQGC-utils/Helpers/get_nanuq_files.py --run ${FC_SHORT}
            if [[ -f "${WORKDIR}/${fc}/SampleSheet.csv" ]]; then
                echo "Launching BCL-convert for run in ${dir}/${fc}/..."
                qsub /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_launcher.sh ${fc}
            else
                echo "ERROR: SampleSheet.csv not found in ${WORKDIR}/${fc}" >&2
                exit 1
            fi
        fi
    else
        echo "Waiting for ${dir}/${fc}/CopyComplete.txt"
    fi
}

for dir in ${WATCHDIRS[@]}; do
    echo "Scanning ${dir}..."
    for FC in $( ls ${dir} ); do
        parts=($(echo ${FC} | tr '_' '\n'))
        if [[ ${#parts[@]} -eq 4 ]]; then
            # Check SampleSheet if run is LowPass.
            # If no SampleSheet is found, then run is not LowPass
            if [[ -f "${dir}/${FC}/SampleSheet.csv" ]]; then
                if grep -q "LowPass" "${dir}/${FC}/SampleSheet.csv"; then
                    echo "PASS: ${dir}/${FC}: Found the word LowPass in SampleSheet"
                else
                    echo "${dir}/${FC}: SampleSheet exists and not for LowPass."
                    launch_run ${dir} ${FC}
                fi
            elif [[ -f "${dir}/${FC}/LowPass*.csv" ]]; then
                echo "PASS: ${dir}/${FC}: Found what looks like a LowPass SampleSheet."
            else
                echo "${dir}/${FC}: Could not find ${dir}/${FC}/SampleSheet.csv."
                    launch_run ${dir} ${FC}
            fi
        fi
        # else: ignore because format of folder name doesn't look like a run
    done
    echo
done
exit

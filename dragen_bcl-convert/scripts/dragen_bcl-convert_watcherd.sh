#!/bin/bash

# Watch for new sequencing runs to launch DRAGEN BCL-Convert and extract stats.
# USAGE (see end of this script for systemd):
#       sudo systemctl start dragen_bcl-convert_watcherd.service
#       nohup /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_watcherd.sh > /dev/null 2>&1 &

BASEDIR="/staging/hiseq_raw"
WORKDIR="/staging2/dragen"
NAPTIME=5
INSTRUMENTS=(A00516 A00977 LH00336 LH00207R)

umask 0002
cd /

while true; do
    # 'CopyComplete.txt' in the run folder marks that sequencing is complete
    # and that all files have been copied to ${BASEDIR}.
    # Launch DRAGEN BCL-Convert if the run folder contains 'CopyComplete.txt', 
    # and if results have not been generated, i.e. output folder exists.
    for instr in ${INSTRUMENTS[@]}; do
        for path_to_run in $( ls -d ${BASEDIR}/${instr}/* ); do
            run=$( echo ${path_to_run} | cut -d / -f 5 )
            if [ -f ${path_to_run}/CopyComplete.txt ]; then
                if [ -d ${WORKDIR}/${run} ]; then
                    # Extract stats from demux
                    if [ -f "${WORKDIR}/${run}/1.fastq/Reports/Demultiplex_Stats.csv" ]; then
                        if [ "${1}" = 'debug' ]; then echo "Extracting demux stats from ${WORKDIR}/${run}/1.fastq/Reports/Demultiplex_Stats.csv"; fi
                        # TODO: python script.py ${run}
                    fi
                else
                    # TODO: Check SampleSheets ? 
                    if [ "${1}" = 'debug' ]; then echo "qsub dragen_bcl-convert_launcher.sh ${run}"; fi
                    # qsub dragen_bcl-convert_launcher.sh ${run}
                    # TODO: logger to file instead of default /var/log/messages ?
                    # touch /var/log/bcl-convert/dragen_bcl-convert_watcherd.log
                    # logger -f /var/log/bcl-convert/dragen_bcl-convert_watcherd.log "Submitted demux job for ${run}"
                    logger "Submitted demux job for ${run}"
                fi
           fi
        done
    done
    if [ "${1}" = 'debug' ]; then echo "Napping for ${NAP}..."; fi
    sleep ${NAPTIME}
done

## Using systemd:
# sudo chmod +x /staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_watcherd.sh
# sudo vi /etc/systemd/system/dragen_bcl-convert_watcherd.service
# Example content of a service file:
#
# [Unit]
# Description=Watch sequencing BCL folders for new runs to demux

# [Service]
# ExecStart=/staging2/soft/CQGC-utils/Analysis.dragen_bcl-convert/scripts/dragen_bcl-convert_watcherd.sh
# Restart=always
# User=root
# StandardOutput=null
# StandardError=null

# [Install]
# WantedBy=multi-user.target
#
# sudo systemctl start myexample.service
# sudo systemctl stop myexample.service
# sudo systemctl enable myexample.service # To start at boot
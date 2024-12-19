#!/bin/bash

# Watch for new sequencing runs to launch DRAGEN BCL-Convert and extract stats.
# USAGE (see end of this script for systemd):
#       sudo systemctl start dragen_bcl-convert_watcherd.service
#       nohup path/to/script.sh > /dev/null 2>&1 &

BASEDIR="/staging/hiseq_raw"
WORKDIR="/staging2/dragen"
INSTRUMENTS=(A00516 A00977 LH00336 LH00207R)
NAP=5

umask 0002
cd /
# touch /var/log/bcl-convert/dragen_bcl-convert_watcherd.log

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
                    echo "Extracting demux stats from ${WORKDIR}/${run}/1.fastq/Reports/DemultiplexStats.csv"
                else
                    # TODO: Check SampleSheets ? 
                    echo "qsub dragen_bcl-convert_launcher.sh ${run}"
                fi
           fi
        done
    done
    echo "Napping for ${NAP}..."
    sleep ${NAP}
done

## Using systemd:
# sudo chmod +x /path/to/script.sh
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
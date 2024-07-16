#!/usr/bin/env python3
"""
Upload FASTQ files to BaseSpace.

USAGE: emg_upload_fastqs.py
       emg_upload_fastqs.py --file samples_list.csv
       emg_upload_fastqs.py --help

Use run and samples information contained in `--file=samples_list.csv` to look 
for FASTQs for uploading to BaseSpace. "samples_list.csv" is normally created 
by the script `emg_collect_dragen_metrics.py` which should have been executed
beforehand. "samples_list.csv" has the following 12 columns:

sample_name,biosample,relation,gender,ep_label,mrn,cohort_type,status,
family_id,birthdate,flowcell_date,flowcell

Tokens to connect to BaseSpace (BSSH) is expected to be found in:
~/.illumina/gapp_conf.json (available at https://github.com/CQGC-Ste-Justine/PrivateDoc/)

Options:

--file="samples_list.csv", list of samples in CSV format. Ex:

    sample_name,biosample,relation,gender,ep_label,mrn,cohort_type,status,family_id,birthdate,flowcell_date,flowcell
    GM241567,27556,PROBAND,FEMALE,CHUSJ,03486257,TRIO,AFF,03486257,2024-04-29,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241601,27560,MTH,FEMALE,CHUSJ,03487612,TRIO,UNF,03486257,1980-10-15,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
    GM241575,27559,FTH,MALE,CHUSJ,03487451,TRIO,UNF,03486257,1978-10-02,2024-07-05,20240705_LH00336_0073_A22MFJFLT3
"""

import os
import argparse
import logging
import subprocess
import pandas as pd
from glob import glob

__version__ = "0.1"

# List of PRAGMatIQ projects on BaseSpace, as of 2024-05-16
# BSSH Project Id required for uploading FASTQs to the right project folder.
#
# bs -c cac1 project list --filter-term "PRAGMatIQ"
# +-----------------+---------+----------------+
# |      Name       |   Id    |   TotalSize    |
# +-----------------+---------+----------------+
# | PRAGMatIQ_CHUSJ | 3703702 | 18289900071562 |
# | PRAGMatIQ_CHUS  | 3703703 | 4840873989012  |
# | PRAGMatIQ_CHUQ  | 4714713 | 2470539400235  |
# | PRAGMatIQ_CUSM  | 5412410 | 2181153257963  |
# +-----------------+---------+----------------+
#
# Hash to link Nanuq samples "ep" to BaseSpace "ProjectId"
#
project_ids = {'CHUSJ': '3703702', 
               'CHUS' : '3703703', 
               'CHUQ' : '4714713', 
               'CUSM' : '5412410'}

def parse_args():
    parser = argparse.ArgumentParser(description="Upload FASTQ files to BaseSpace. for a given Run.")
    parser.add_argument('--file', '-f', default="samples_list.csv", help="Get samples from file. Default='samples_list.csv'.")
    parser.add_argument('--logging-level', '-l', dest='level', default='info',
                        help="Logging level (str), can be 'debug', 'info', 'warning'. Default='info'")
    return(parser.parse_args())


def configure_logging(level):
    """
    Set logging level, based on the level names of the `logging` module.
    - level (str): 'debug', 'info' or 'warning'
    """
    if level == 'debug':
        level_name = logging.DEBUG
    elif level == 'info':
        level_name = logging.INFO
    else:
        level_name = logging.WARNING
    logging.basicConfig(level=level_name, 
                        format='[%(asctime)s] %(levelname)s: %(message)s', 
                        datefmt='%Y-%m-%d@%H:%M:%S')


def main(args):
    """
    Iterate through information in "samples_list.csv" to build `bs` command for
    uploading FASTQs to BaseSpace.
    """
    df = pd.read_csv(args.file)

    sites = df['Site'].unique()
    logging.info(f"Uploading {len(df)} samples for {len(sites)} site:")
    for site in sites:
        logging.info(f"{site} => {len(df[df['Site'] == site])}")

    for row in df.itertuples():
        logging.info(f"List FASTQs for biosample={row.biosample} to upload to BBSH folder PRGAMatIQ_{row.ep_label}")
        fastqdir = f"/staging/hiseq_raw/{row.flowcell.split('_')[1]}/{row.flowcell}/Analysis/1/Data/DragenGermline/fastq"
        os.chdir(fastqdir)
        fastqs = glob(f"{row.biosample}_*.fastq.gz")
        results = subprocess.run((['bs', '-c', 'cac1', 'dataset', 'upload', 
                                    '--no-progress-bars', 
                                    '--project', f"{project_ids[row.ep_label]}", 
                                    '--biosample-name', f"{row.biosample}"] + fastqs), 
                                    capture_output=True, text=True)
        logging.info(f"Upload to BSSH complete for {row.biosample} (STDOUT):\n{results.stdout}")
        logging.debug(f"stdargs:\n{results.stderr}")
        logging.debug(f"stderr:\n{results.stderr}")
        
        # ```bash
        # ep="CHUSJ"
        # project=$( bs -c cac1 projects list --terse --filter-term "PRAGMatIQ_${ep}$" )
        # for sample in 27556 27560 27559 27555 27557 27558; do
        #     fastqs=$( ls ${sample}_*_001.fastq.gz )
        #     bs -c cac1 dataset upload --project ${project} --biosample-name ${sample} ${fastqs}
        # done
        # ```


def tests(args):
    """
    Write some quick and dirty tests
    """
    print(args)
    return(0)


if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)
    #tests(args)

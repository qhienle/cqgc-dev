#!/usr/bin/env python3
"""
Upload FASTQ files to BaseSpace.

USAGE: emg_upload_fastqs.py
       emg_upload_fastqs.py --file samples_list.csv
       emg_upload_fastqs.py --help

Use run and samples information contained in `--file=samples_list.csv` to look 
for FASTQs for uploading to BaseSpace. "samples_list.csv" is normally created 
by the script `list_run_samples.py` which should have been executed beforehand.
File "samples_list.csv" has the following 11 columns:

sample_name,biosample,relation,gender,ep_label,mrn,status,family_id,birthdate,project,flowcell

Tokens to connect to BaseSpace (BSSH) is expected to be found in:
~/.illumina/gapp_conf.json (available at https://github.com/CQGC-Ste-Justine/PrivateDoc/)

Options:

--file="samples_list.csv", list of samples in CSV format. Ex:

sample_name,biosample,relation,gender,ep_label,mrn,status,family_id,birthdate,project,flowcell
Q1K_HSJ_1525-1130_S2,36217,SIS,FEMALE,CHUSJ,0,AFF,1525-1130,2019-12-12,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
Q1K_HSJ_1525-1130_P,36218,PROBAND,MALE,CHUSJ,0,AFF,1525-1130,2013-06-23,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
Q1K_HSJ_1525-1130_M1,36215,MTH,FEMALE,CHUSJ,0,UNF,1525-1130,1987-04-12,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
Q1K_HSJ_1525-1130_S1,36216,BRO,MALE,CHUSJ,0,UNF,1525-1130,2017-07-05,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
"""

import os
import sys
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
# Hash to link Nanuq samples "ep" to BaseSpace "ProjectId". This is the project
# folder where FASTQ files are stored. Projects Q1K and AOH do not distinguish
# between EP_Labels (Établissement Public) like PRAGMatIQ.
#
project_ids = {'CHUSJ': '3703702', 
               'CHUS' : '3703703', 
               'CHUQ' : '4714713', 
               'CUSM' : '5412410',
               'Q1K'  : '6197214',
               'AOH'  : '6050046'}

def parse_args():
    parser = argparse.ArgumentParser(description="Upload FASTQs to BaseSpace for samples listed in --file(=samples_list.csv).")
    parser.add_argument('--file',     '-f', default="samples_list.csv", help="Get samples from file. Default='samples_list.csv'.")
    parser.add_argument('--data-dir', '-d', help="Get FASTQs from --data-dir. Default='fastq folder'.")
    parser.add_argument('--project',  '-p', default='prag', help="Project: 'prag', 'eval', 'q1k', 'aoh'. Default='prag'")
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

    logging.info(f"Uploading {len(df)} samples to BaseSpace :")

    # Set Project ID for PRAGMatIQ EP Labels. AOH and Q1K store all FASTQs of a
    # project in the same folder. PRAGMatIQ (incl. eval) separates between EP.
    #
    if args.project == 'prag' or args.project == 'eval':
        logging.info(f"Not updating 'ep_label' in 'samples_list.csv' for project PRAGMatIQ")
    elif args.project == 'q1k':
        df['ep_label'] = 'Q1K'
    elif args.project == 'aoh':
        df['ep_label'] = 'AOH'
    else:
        logging.error(f"Bad project name {args.project}. Must be: prag, eval, q1k or aoh")
        sys.exit()
    for ep in df['ep_label'].unique(): logging.info(f"{ep} => {len(df[df['ep_label'] == ep])}")

    # List FASTQ files for each sample and upload to BaseSpace
    #
    for row in df.itertuples():
        logging.info(f"List FASTQs for biosample={row.biosample} to upload to BBSH folder PRGAMatIQ_{row.ep_label}")
        if args.data_dir is not None:
            fastqdir = args.data_dir
        else:
            #fastqdir = f"/staging2/dragen/{row.flowcell}/1.fastq"
            fastqdir = f"/mnt/vs_nas_chusj/CQGC_PROD/fastqs/{row.flowcell}/1.fastq"
        os.chdir(fastqdir)

        # glob() does not create ordered list of files, so we sort()for `bs`
        # `bs` checks for R1/R2 pairs and panics if both files are not listed
        # consecutively.
        #
        fastqs = glob(f"{row.biosample}_*.fastq.gz")
        fastqs.sort() 
        results = subprocess.run((['bs', '-c', 'cac1', 'dataset', 'upload', 
                                    '--no-progress-bars', 
                                    '--project', f"{project_ids[row.ep_label]}", 
                                    '--biosample-name', f"{row.biosample}"] + fastqs), 
                                    capture_output=True, text=True)
        if results.stderr != '':
            logging.warning(f"ERROR while subprocess.run():\n{results.stderr}")
            logging.warning(f"args:\n{results.args}")
        else:
            logging.info(f"Upload to BSSH complete for {row.biosample} (STDOUT):\n{results.stdout}")
        
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

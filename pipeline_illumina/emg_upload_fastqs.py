#!/usr/bin/env python3
"""
Upload FASTQ files to BaseSpace.

USAGE: emg_upload_fastqs.py
       emg_upload_fastqs.py --file samples_list.csv
       emg_upload_fastqs.py --help

Options:

--file, list of samples in CSV format. Default

Tokens to connect with Phenotips and BaseSpace (BSSH) are expected to be found 
in ~/.illumina/gapp_conf.json (available at https://github.com/CQGC-Ste-Justine/PrivateDoc/)
"""

import os, sys
import argparse
import logging
import json
import re
import subprocess
import pandas as pd

# Set source path to CQGC-utils so that we can use relative imports
#
src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(src_path)
from lib.nanuq import Nanuq
from lib.gapp import BSSH

nq = Nanuq()

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
# This hash to link "ep" fields for samples in Nanuq to BSSH ProjectId
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
    """
    try: 
        samples_list = pd.read_csv('samples_list.csv')
    except FileNotFoundError as err:
        logging.error(err)
    else:
        pass
    
    fastqdir = f"/staging/hiseq_raw/{fc_short}/{args.run}/Analysis/1"
    try:
        os.mkdir(workdir)
        logging.info(f"Created work directory '{workdir}'")
    except:
        logging.warning(f"Could not create {workdir}")

    # 1. Get a list of samples on this run to construct the cases.
    # TODO: Add experiment name as an alternative identifier for Nanuq API?
    #
    df = nanuq_samples_to_df(args.run)
    df['fc_date'] = fc_date
    logging.info(f"Add column for flowcell date {fc_date}")
    df = df.sort_values(by=['Family Id', 'relation'], ascending=[True, False])
    df.to_csv(f"{workdir}{os.sep}samples_list.csv", index=None)
    logging.info(f"Saved list of samples to file: {workdir}{os.sep}samples_list.csv")
    
    for biosample in df['biosample']: 
        logging.info(f"Uploading FASTQs for {biosample}")


def tests(args):
    return(1)

if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)
    #tests(args)


"""
Example of Nanuq SampleNames:
```pyhton
nq = Nanuq()
samplenames = nq.get_samplenames("A00516_0445")
samplenames.text
```
##2023-08-09
##Centre for Pediatric Clinical Genomics
##Flow Cell: H7N33DSX7
##Principal Investigator: Dr Mehdi Yeganeh, CHUS Service de génétique médicale, Dr Jean-François Soucy
##Nanuq References: A00516_0445
##Content: Internal_Sample_ID -> Client_Sample_Name Conversion grid
##-------------------------------------------
##Internal_Sample_ID	Client_Sample_Name
##-------------------------------------------
22256	3042652455
22257	3042642360
22258	3042645886
22262	3052768938
22263	3052768942
22264	3052768961
22265	3052768951
22282	23-05982-T1
22283	23-06383-T1
22284	23-06384-T1
22285	GM231615
22286	GM231624
22287	GM231626
22288	GM231627
22290	GM231632
22293	GM231651
...
##The "Internal_Sample_ID" is the internal identifier assigned by the Center for this sample.
##Use of an "Internal_Sample_ID" ensures the traceability of this sample throughout the sequencing
##process, but also the anonymization of the information generated and transferred in the
##frame of this project.


Example of Nanuq sample: `json.loads(nq.get_sample(22283))`
[{'ldmSampleId': '23-06383-T1',
  'ldm': 'LDM-CHUS',
  'patient': {'designFamily': 'TRIO',
   'familyId': '23-05982-T1',
   'familyMember': 'FTH',
   'firstName': 'Andre-Philippe',
   'lastName': 'Belley',
   'ramq': 'BELA93092213',
   'sex': 'MALE',
   'mrn': '00000000',
   'ep': 'CHUS',
   'birthDate': '22/09/1993',
   'fetus': False,
   'status': 'UNK'},
  'ldmServiceRequestId': '23-06383-T1',
  'labAliquotId': '22283',
  'panelCode': 'PRAGMATIQ',
  'specimenType': 'NBL',
  'sampleType': 'DNA',
  'ldmSpecimenId': '23-06383-T1'}]
"""
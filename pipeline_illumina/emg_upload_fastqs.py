#!/usr/bin/env python3
"""
Upload FASTQ files to BaseSpace.

USAGE: emg_upload_fastqs.py <RUN>
       emg_upload_fastqs.py 20240510_LH00336_0043_A22K5KMLT3
       emg_upload_fastqs.py --help

List of samples can be provided using --file, instead of fetching from Nanuq.
This file can either be the "SampleNames.txt" downloaded from Nanuq, or a one-
column listing of CQGC IDs.

Nanuq username and password have be saved in a file named '~/.nanuq', like so:
`echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`
Replace USERNAME and PASSWORD with actual values.

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
bssh = BSSH()

__version__ = "0.1"


def parse_args():
    parser = argparse.ArgumentParser(description="Upload FASTQ files to BaseSpace. for a given Run.")
    parser.add_argument('run', help="Run ID, ex: '20240510_LH00336_0043_A22K5KMLT3'")
    parser.add_argument('--file', '-f', help="Get samples from `--file` instead of fetching <Run> from Nanuq.")
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


def list_samples(file=None):
    """
    Return a list of CQGC ID for samples from `file`. If file is None, get list
    from Nanuq directly.
    - file (str): Nanuq's "SampleNames.txt" or a one-column list of CQGC IDs.
    """
    samples = []
    if file is None:
        samplenames = nq.get_samplenames(args.run)
        if not samplenames.text.startswith("##20"):
            sys.exit(logging.error(f"Unexpected content for SampleNames. Please verify Nanuq's reponse:\n{samplenames.text}"))
        else:
            logging.info("Retrieved samples conversion table from Nanuq")
            lines = samplenames.text.splitlines()
    else:
        logging.info(f"Using list of samples from file {args.file} instead of Nanuq")
        with open(file, 'r') as fh:
            lines = fh.readlines()

    for line in lines:
        if not line.startswith('#'):
            samples.append(line)
    return samples


def main(args):
    """
    """
    date, inst, args.run.split('_')
    print(f"# Logging run {args.run}")
    
    workdir = f"{os.getcwd()}{os.sep}{args.run}"
    try:
        os.mkdir(workdir)
        logging.info(f"Created work directory '{workdir}'")
    except:
        logging.warning(f"Could not create {workdir}")

    # 1. Get a list of samples on this run to construct the cases.
    # TODO: Add experiment name as an alternative identifier for Nanuq API?
    #
    samplenames = list_samples(args.file)
    for line in samplenames:
        cqgc, sample = line.split("\t")
        print(cqgc)


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
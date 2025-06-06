#!/usr/bin/env python3
"""
Create cases in Emedgene from a file listing individual and family informations 
for samples (default='samples_list.csv').

USAGE: emg_create_cases.py
       emg_create_cases.py -f samples_list.cav
       emg_create_cases.py --help

Tokens to connect with Phenotips and BaseSpace (BSSH) are expected to be found 
in ~/.illumina/gapp_conf.json (available at https://github.com/CQGC-Ste-Justine/PrivateDoc/)
"""

import os, sys
import argparse
import logging
import json
import requests
import pandas as pd

# Set source path to CQGC-utils so that we can use relative imports
#
src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(src_path)
from lib.gapp import BSSH
from lib.gapp import Emedgene
from lib.nanuq import Nanuq
bssh = BSSH()

__version__ = "0.2"

def parse_args():
    parser = argparse.ArgumentParser(description="Create cases in Emedgene from list of samples with information for cases")
    parser.add_argument('--file',    '-f', nargs='?', default='samples_list.csv', help="List of samples with Case information. Default='samples_list.csv'")
    parser.add_argument('--project', '-p', default='prag', help="Project: 'prag', 'eval', 'q1k', 'aoh'. Default='prag'")
    parser.add_argument('--logging-level', '-l', dest='level', default='info',
                        help="Logging level (str), can be 'debug', 'info', 'warning'. Default='info'")
    return parser.parse_args()


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


def get_hpos(pid, db='D:\\HSJ\\Workspace\\cqgc-dev\\lib\\pheno_json-extracted_hpos.json'):
    """
    Lookup Phenotips ID (PID) and HPO identifiers from `db`
    - db: [str] Path to JSON file containing PID to HPO information
    - Returns:
    """
    with open(db, "r") as fh:
        pids = json.load(fh)
    try:
        pids[pid]
    except KeyError:
        return('')
    else:
        hpos = []
        for hpo in pids[pid]:
            hpos.append(hpo['id'])
        return(';'.join(hpos))


def get_birthdate(biosample):
    """
    Get date of birth from Nanuq for `biosample`
    """
    nq = Nanuq()
    data = json.loads(nq.get_sample(biosample))
    return(data[0]['patient']['birthDate'])


def submit_emg(case_json):
    """
    Submit to Emedgene a Case JSON. 
    - case_json: [str] path/to/json
    - Return:    [obj] Response object from a `requests()` POST.
    """
    emg = Emedgene()
    with open(case_json, 'r') as fh:
        url = f"{emg.eval_server}/api/cases/v2/cases/"
        payload = json.load(fh)
        resp = requests.post(url, headers={'Authorization': emg.authenticate()}, json=payload)
    return resp.json()


def main():
    """
    1. Get individual and family information from a samples' list (csv)
    2. Add BaseSpace FASTQ file paths for each sample
    3. Get the corresponding HPO Identifiers and add HPO terms
    4. Add storage provider and label ID's based on project
    5. Write Emedgene case JSON, including paths to VCFs and HPO terms;
    6. Create case on Emedgene
    """
    args = parse_args()
    configure_logging(args.level)
    workdir = os.getcwd()
    os.chdir(workdir)

    # 1. Get individual and family information from a samples' list (csv)
    try:
        df_samples = pd.read_csv(args.file)
    except FileNotFoundError as e:
        logging.error(e)
    except Exception as e:
        logging.error(f"Could not load list of samples from '{args.file}' because {e}")
    else:
        df_samples['mrn'] = df_samples['mrn'].astype(str)
        workdir = os.path.dirname(os.path.abspath(args.file))
        os.chdir(workdir)
        logging.info(f"Log {__file__} for run {','.join(df_samples['flowcell'].unique())}")
        logging.debug(f"Run arguments: {args}")


    # 2. Add BaseSpace FASTQ file paths for each sample
    # 3. Get the corresponding HPO Identifiers and add HPO terms
    # 4. Add storage provider and label ID's based on project
    # 5. Write Emedgene case JSON, including paths to VCFs and HPO terms;
    # 6. Create case on Emedgene


if __name__ == '__main__':
    sys.exit(main())

"""
## EXAMPLE: A family in 'samples_list.csv'

sample_name,biosample,relation,gender,ep_label,mrn,status,family_id,birthdate,project,flowcell
Q1K_HSJ_1525-1130_S2,36217,SIS,FEMALE,CHUSJ,0,AFF,1525-1130,2019-12-12,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
Q1K_HSJ_1525-1130_P,36218,PROBAND,MALE,CHUSJ,0,AFF,1525-1130,2013-06-23,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
Q1K_HSJ_1525-1130_M1,36215,MTH,FEMALE,CHUSJ,0,UNF,1525-1130,1987-04-12,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4
Q1K_HSJ_1525-1130_S1,36216,BRO,MALE,CHUSJ,0,UNF,1525-1130,2017-07-05,Q1K_CHUSJ,20250523_LH00336_0218_B22TNY2LT4

## Alternative workflow:

0. Setup run and Quality Check
1. Get list of samples from Nanuq or from file in order to construct cases.
    For each sample:
    1.1 Get family information from Nanuq
    1.2. Add the Phenotips ID (PID) and the corresponding HPO Identifiers;
2. Run joint joint variant calling for each case;
3. Upload the VCF output from joint genotyping to EMedgene AWS S3 bucket;
4. Write Emedgene case JSON, including paths to VCFs and HPO terms;
5. Create case on Emedgene

Example for joint-genotyping:
/staging2/soft/CQGC-utils/Analysis.pipeline_exome/pipeline.exomes.DRAGEN.joint-genotyping.sh
"""

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
import re
import pandas as pd

# Set source path to CQGC-utils so that we can use relative imports
#
src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(src_path)
from lib.gapp import BSSH
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


def list_samples_from_samplenames(content):
    """
    Get list of samples from the content of Nanuq's "SampleNames.txt".
    A one-column list of CQGC identifiers will also work.
    - content: [str] File content of SampleNames.txt, as a single string.
    - Returns: [list] of samples from Nanuq (CQGC IDs, referred to by 
               Illumina as 'biosamples')
    """
    biosamples = []
    for line in content.splitlines():
        logging.debug(f"Parsing SampleNames...")
        if not line.startswith("#"):
            try:
                cqgc_id, sample_name = line.split("\t")
            except ValueError as err:
                logging.debug(f"While parsing SampleNames content: {err}.")
            else:
                biosamples.append(cqgc_id)
    return biosamples


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


def main():
    """
    0. Setup run and Quality Check
    1. Get list of samples from Nanuq or from file in order to construct cases.
       For each sample:
        1.1 Get family information from Nanuq
        1.2. Add the Phenotips ID (PID) and the corresponding HPO Identifiers;
    2. Run joint joint variant calling for each case;
    3. Upload the VCF output from joint genotyping to EMedgene AWS S3 bucket;
    4. Write Emedgene case JSON, including paths to VCFs and HPO terms;
    5. Create case on Emedgene
    """

    # 0. Setup for run
    #
    args = parse_args()
    configure_logging(args.level)
    workdir = os.getcwd()
    os.chdir(workdir)

    nq = Nanuq()

    # 1. Get list of cases and family information from Nanuq or SampleNames.txt
    #
    # biosamples = list_samples_from_samplenames(args.run, file=args.file)
    # logging.info(biosamples)
    if args.run:
        logging.info(f"Listing samples from Nanuq for run {args.run}.")
        biosamples = list_samples_from_samplenames(nq.get_samplenames(args.run).text)
    elif args.file:
        logging.info(f"Listing samples from file {args.file}.")
        with open(args.file, 'r') as fh:
            biosamples = list_samples_from_samplenames(fh.read())
    logging.info(f"List of {len(biosamples)} samples to process: {biosamples}")

    # 1.1 Get family information from Nanuq. Parse and store as a list of data 
    # structures in "cases", that can be loaded into a pandas DataFrame for
    # easier data wrangling.
    #

    cases = []
    for sample in biosamples:
        try:
            data = json.loads(nq.get_sample(sample))
        except Exception as e:
            logging.warning(f"JSONDecodeError {e} could not decode biosample {sample}")
            continue

        logging.info(f"Got information for biosample {sample}")
        if len(data) != 1:
            logging.debug(f"Number of samples retrieved from Nanuq is not 1.\n{data}")
        try:
            data[0]["patient"]["mrn"]
        except Exception as err:
            logging.warning(f"Could not find MRN for patient {sample}: {err}")
            data[0]["patient"]["mrn"] = '0000000'
        else:
            pass
        finally:
            sample_infos = [
                data[0]["ldmSampleId"],
                data[0]["labAliquotId"],
                data[0]["patient"]["familyMember"],
                data[0]["patient"]["sex"],
                data[0]["patient"]["ep"],
                data[0]["patient"]["mrn"],
                data[0]["patient"]["designFamily"],
                data[0]["patient"]["birthDate"],
                data[0]["patient"]["status"],
                data[0]["patient"].get("familyId", "-")
            ]

        # 1.2 Add Phenotips ID (`pid`) and patients' HPO identifiers for
        # the proband. Lookup this information in Phenotips, using EP+MRN
        # Ex: CHUSJ123456
        #
        if data[0]["patient"]["familyMember"] == 'PROBAND':
            #pid, labels_str, ids_str = add_hpos(data[0]["patient"]["ep"], data[0]["patient"]["mrn"])
            logging.info(f"Got HPO terms from Phenotips for PID {pid}")
        else:
            pid, labels_str, ids_str = ('', '', '')
            logging.debug(f'Not retrieving PID for {sample} ({data[0]["patient"]["familyMember"]})')
        sample_infos.append(pid)
        sample_infos.append(labels_str)
        sample_infos.append(ids_str)
        logging.debug(f"PID: {pid}; HPO ID: {ids_str}; Labels: {labels_str}\n")

        cases.append(sample_infos)

    # For each case, prepare the JSOn payload by replacing the template's info
    # with samples data.
    #
    df = pd.DataFrame(cases)
    df.columns = ['sample_name', 'biosample', 'relation', 'gender', 'label', 
                  'mrn', 'cohort_type', 'date_of_birth(YYYY-MM-DD)', 'status',
                  'Family Id', 'pid', 'phenotypes', 'hpos', 'filenames']
    df = df.sort_values(by=['Family Id', 'relation'], ascending=[True, False])
    print(df)

if __name__ == '__main__':
    main()

"""
## EXAMPLE: A trio in 'samples_list.csv'

sample_name,biosample,relation,gender,ep_label,mrn,status,family_id,birthdate,flowcell_date,flowcell
DM250647,33641,PROBAND,FEMALE,CHUSJ,03514485,AFF,03514485,2024-04-08,2025-02-12,20250212_LH00336_0170_B22V7KMLT3
DM250631,33643,MTH,FEMALE,CHUSJ,03381944,UNF,03514485,1992-10-12,2025-02-12,20250212_LH00336_0170_B22V7KMLT3
DM250630,33642,FTH,MALE,CHUSJ,01483718,UNF,03514485,1984-01-27,2025-02-12,20250212_LH00336_0170_B22V7KMLT3


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

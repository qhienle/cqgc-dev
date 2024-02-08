#!/usr/bin/env python3
"""
Create cases in Emedgene for <run> (ex.: "A00516_420").

USAGE: emg_create_cases.py [-f] <run>
       emg_create_cases.py A00516_420
       emg_create_cases.py A00516_420 --file SampleNames.txt
       emg_create_cases.py --help

Create cases in Emedgene for samples listed for <run> (_e.g._ LH00336_0009).
Sample information is downloaded from Nanuq based on the CQGC ID. 
List of samples can be provided using --file. This file can either be the 
"SampleNames.txt" downloaded from Nanuq, or a one-column listing of CQGC IDs.

Nanuq username and password have be saved in a file named '~/.nanuq', like so:
`echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`
Replace USERNAME and PASSWORD with actual values.

Tokens to connect with Phenotips and BaseSpace (BSSH) are expected to be found 
in ~/.illumina/gapp_conf.json (available at https://github.com/CQGC-Ste-Justine/PrivateDoc/)

Example for joint-genotyping:
/staging2/soft/CQGC-utils/Analysis.pipeline_exome/pipeline.exomes.DRAGEN.joint-genotyping.sh
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
    parser = argparse.ArgumentParser(description='Create cases in Emedgene for <run> (ex.: "A00516_420").')
    parser.add_argument('run', help="FC_SHORT Run ID, ex: 'A00516_339'")
    parser.add_argument('--file', '-f', help="List of samples", default='SampleNames.txt')
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


def list_samples_from_samplenames(content):
    """
    Get list of samples from the content of Nanuq's "SampleNames.txt".
    A one-column list of CQGC identifiers will also work.
    - content: [str] File content of SampleNames.txt, as a single string.
    - Returns: [list] of samples from Nanuq (CQGC IDs, referred to by 
               Illumina as 'biosamples')
    """
    biosamples = []
    for line in content:
        logging.debug(f"Parsing SampleNames...")
        if line.startswith("#"):
            if line.startswith("##20"):
                # Is this still useful?
                fc_date = re.match(r'##(\d{4}-\d{2}-\d{2})', line).group(1)
        else:
            try:
                cqgc_id, sample_name = line.split("\t")
            except ValueError as err:
                logging.debug(f"While parsing {file}: {err}.")
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
    1. Get list of samples and family information from Nanuq or SampleNames.txt
       in order to construct cases;
    2. Run joint joint variant calling for each case;
    3. Upload the VCF output from joint genotyping to EMedgene AWS S3 bucket;
    4. Write Emedgene case JSON;
        4.1. Include paths to VCFs and HPO terms;
        4.2. Add the Phenotips ID (PID) and the corresponding HPO Identifiers;
    5. Create case on Emedgene
    """

    args = parse_args()
    configure_logging(args.level)
    workdir = os.getcwd()
    os.chdir(workdir)

    # 1. Get list of cases and family information from Nanuq or SampleNames.txt
    #
    biosamples = list_samples_from_samplenames(args.run, file=args.file)
    logging.info(biosamples)


if __name__ == '__main__':
    main()

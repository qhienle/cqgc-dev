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


def list_samples_from_samplenames(run, file=None):
    """
    Get list of samples, either directly from Nanuq (default) or from the Nanuq
    file, "SampleNames.txt"
    - db: [str] Path to JSON file containing PID to HPO information
    - Returns: List of samples from Nanuq (CQGC IDs, referred to by 
               Illumina as 'biosamples')
    """
    biosamples = []
    nq = Nanuq()
    if file is not None:
        with open(file, 'r') as fh:
            content = fh.read().splitlines()
        logging.info("Listing samples... Retrieved samples conversion table from file {file}.")
    else:
        content = nq.get_samplenames(run).text.splitlines()
        logging.info("Listing samples... Retrieved samples conversion table from Nanuq.")

    for line in content:
        logging.debug(f"Parsing SampleNames...")
        if line.startswith("#"):
            if line.startswith("##20"):
                fc_date = re.match(r'##(\d{4}-\d{2}-\d{2})', line).group(1)
        else:
            cqgc_id, sample_name = line.split("\t")
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


def df_to_manifest(df):
    """
    From data in df, generate a manifest file for batch upload to Emedgene 
    using the UI. For specifications of the manifest, see:
    "https://help.emedgene.com/en/articles/7231644-csv-format-requirements"
    - `df`: A Pandas DataFrame
    - Returns: File 'emg_batch_manifest.csv' in current folder
    """
    df_manifest = pd.DataFrame({
        'Family Id': df['family'],
        'Case Type': 'Whole Genome',
        'Files Names': df['filenames'],
        'Sample Type': 'FASTQ',
        'BioSample Name': df['sample_name'],
        'Visualization Files': '',
        'Storage Provider Id': 10123,
        'Default Project': '',
        'Execute_now':  'False',
        'Relation': df['relation'],
        'Gender': df['gender'],
        'Phenotypes': 'Healthy',
        'Phenotypes Id': df['hpos'],
        'Date Of Birth': pd.to_datetime(df['birthdate'], format='%d/%m/%Y'),
        'Boost Genes': '',
        'Gene List Id': '',
        'Kit Id': '',
        'Selected Preset': 'Default',
        'Label Id': df['label'],
        'Clinical Notes': df['pid'],
        'Due Date': '',
        'Opt In': ''
    })
    # With the "Files Names"="auto" option, BSSH users can automatically locate
    # FASTQ files based on the BioSample Name and Default Project provided.
    # Unfortunately, this would mean that cases woul bear the lab's CQGC_ID.
    #
    df_manifest.loc[df_manifest['Relation'] == 'PROBAND', 'Phenotypes'] = ''
    df_manifest['Default Project'] = 'Rapidomics'

    df_manifest['Relation'].replace('PROBAND', 'proband', inplace=True)
    df_manifest['Relation'].replace('MTH', 'mother', inplace=True)
    df_manifest['Relation'].replace('FTH', 'father', inplace=True)
    df_manifest['Relation'].replace('BRO', 'sibling', inplace=True)
    df_manifest['Relation'].replace('SIB', 'sibling', inplace=True) # TODO: Verify 'SIB', or SIS?

    df_manifest['Gender'].replace('FEMALE', 'F', inplace=True)
    df_manifest['Gender'].replace('MALE', 'M', inplace=True)
    df_manifest['Gender'].replace('', 'U', inplace=True) 

    with open('emg_batch_manifest.csv', 'w') as fh:
        fh.write('[Data],,,,,,,,,,,,,,,,,,,,,\n')
        fh.write(df_manifest.to_csv(index=None, lineterminator='\n'))


def main():
    """
    1. Get list of cases and family information from Nanuq or SampleNames.txt;
    2. Add the Phenotips ID (PID) and the corresponding HPO Identifiers;
    3. Connect to BaseSpace and re-construct the path to the FASTQ files;
    4. Convert DataFrame into a CSV file (manifest) for EMG batch upload.
    """

    args = parse_args()
    configure_logging(args.level)
    workdir = os.getcwd()
    os.chdir(workdir)

    # 1. Get list of cases and family information from Nanuq or SampleNames.txt
    #


if __name__ == '__main__':
    main()

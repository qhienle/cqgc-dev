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
    db_xlsx = args.file
    workdir = 'D:\HSJ\Projects\PRAGMatIQ\TSS Negatives AOH'
    os.chdir(workdir)

    # 1. Get list of cases and family information as a pandas dataframe.
    #
    logging.debug(f"Reading file {db_xlsx}")
    df_rapidomics0 = pd.read_excel(db_xlsx, sheet_name='Rapidomics_Samples', usecols='A:F,G,I')
    df_rapidomics0.columns = ['sample_name', 'biosample', 'gender', 'pid', 'relation', 'design', 'family', 'run']
    df_unsolved0   = pd.read_excel(db_xlsx, sheet_name='Unsolved', usecols='F')
    df_unsolved0.columns = ['sample_name']

    df = pd.merge(df_unsolved0, df_rapidomics0, how='inner', on=['sample_name'])

    # 2. Add the Phenotips ID (PID) and the corresponding HPO Identifiers
    # 
    df['hpos'] = df['pid'].apply(get_hpos)

    # 3. Connect to BaseSpace and re-construct the path to the FASTQ files
    #
    df['filenames'] = ''
    print(f"Path to 13203 {bssh.get_sequenced_files('13203')}")    

    # 4. Convert DataFrame into a CSV file (manifest) for EMG batch upload
    #
    df['label'] = 15 # LabelID=15 for CHUSJ, on eval
    df['birthdate'] = df['biosample'].apply(get_birthdate)
    df = df.sort_values(by=['family', 'relation'], ascending=[True, False])
    print(df)
    df_to_manifest(df)
    logging.info("Wrote manifest file `emg_batch_manifest.csv` for batch upload to Emedgene.")


if __name__ == '__main__':
    main()

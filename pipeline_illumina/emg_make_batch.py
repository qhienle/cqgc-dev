#!/usr/bin/env python3
"""
Make a batch file for Case creation in Emedgene from list of samples.

USAGE: emg_make_batch_from_nanuq.py --file samples_list.csv
       emg_make_batch_from_nanuq.py --help

List of samples can either be the "SampleNames.txt" downloaded from Nanuq, 
or a one-column listing of CQGC IDs.

Credentials to access Emedgene, Phenotips and REDCap should be stored in a 
configuration file in JSON format. The default `gapp_conf.json` must contain:

{
    "instance"         : "cac1.trusight.illumina.com",
    "X-ILMN-Domain"    : "chusj",
    "X-ILMN-Workgroup" : "<replace_with_your_hash_key>",
    "X-Auth-Token"     : "APIKey <replace_with_your_hash_key>",
    "testDefinitionId" : "<replace_with_your_hash_key>",
    "bs_apiServer"     : "https://api.cac1.sh.basespace.illumina.com",
    "bs_accessToken"   : "<replace_with_your_hash_key>",
    "X-Gene42-Server"  : "https://chusj.phenotips.com",
    "X-Gene42-Auth"    : "Basic replace_with_your_hash_key",
    "X-Gene42-Secret"  : "<replace_with_your_hash_key>",
    "REDCap-Server"    : "https://tacc-redcap.bic.mni.mcgill.ca/api/",
    "REDCap-Token"     : "<replace_with_your_hash_key>"
}

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
from lib.gapp import Phenotips
from lib.gapp import REDCap
from lib.gapp import BSSH

__version__ = "0.1"


def parse_args():
    parser = argparse.ArgumentParser(description="Make Emedgene batch file for case creation from samples_list.")
    parser.add_argument('--file',    '-f', default='samples_list.csv', help="List of samples with Case information")
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


def add_fastqs(biosample):
    """
    Add BSSH paths to fastq files for samples listed in df
    - biosample: [str] Name of biosample (usually CQGC LabID, ex: 27692)
    - Returns  : [str] df appened with column 'filenames'
    """
    bssh = BSSH()
    try:
        fastqs = bssh.get_sequenced_files(biosample)
    except Exception as err:
        logging.info(f"Could not retrieve FASTQs paths for {biosample}: {err}")
        fastqs = []
    else:
        filenames = ';'.join(fastqs)
    finally:
        return filenames


def add_hpos_phenotips(ep, mrn):
    """
    Lookup Phenotips ID (PID) and HPO identifiers
    - ep     : [str] Etablissement Public. Ex: CHUSJ
    - mrn    : [str] Medical Record Number. Ex: 123456
    - Returns: [str] Semi-column-spearated list of hpo identifiers
    """
    pho = Phenotips()
    pid = ''
    hpo_ids    = []
    hpo_labels = []

    # Call API with labeled external ID (eid) to retrieve PID and then get the 
    # associated HPO terms. In Phenotips, MRN is prepended with the Site's EP 
    # initials (_e.g._ "CHUS1626861"). There does not seem to be a standard 
    # format for MRN identifiers. Here are some of the format detected:
    #
    # Phenotips     | Nanuq            | Notes
    # --------------|------------------|------------------------------------
    # CHUSJ3421069  | 03421069         | For CHUSJ, Nanuq adds a leading "0"
    # CHUSJX3627954 | X3627954         | Not numerical, starts with 'X'
    # CHUS1628699   | 1628699          | 7 digits, no leading '0'
    # CHUS347990    | 347990           | 6 and no leading '0' added by Nanuq
    # 1633799       | 1633799          | Not prefixed with EP initials
    # CHUQ1644460   | CHUL1644460      | CHUQ is stored as CHUL in Nanuq
    # CHUQ1753303   | CHUQ1753303 CHUL | Suffix CHUL added in Nanuq
    #
    # Fix malformed entries.
    #
    if ep == 'CHUSJ':
        mrn = mrn.lstrip('0')
    elif ep == 'CHUS':
        pass
    elif ep == 'CHUL':
        mrn.replace('L', 'Q')
    elif ep == 'MUHC':
        ep = 'CUSM'
    elif mrn.startswith('MCH_'):
        mrn = mrn.replace('MCH_', '')
    elif mrn.endswith(' CHUL'):
        mrn = mrn.replace(' CHUL', '')
    ep_mrn = f"{ep}{mrn}"
    patient  = pho.get_patient_by_mrn(ep_mrn)

    # EP+MRN is a convention, not a constraint enforced in Phenotips DB
    # Users sometimes don't follow the rule and provide only MRN
    # Use RAMQ (not always available) as a last resort? 
    #
    if patient is not None:
        pid = patient['id']
    else:
        logging.warning(f"Could not get PID using EP+MRN: {ep_mrn}. Trying with MRN: {mrn}...")
        if ep == 'CHUSJ':
            patient = pho.get_patient_by_mrn(mrn.lstrip('0')) # Why is MRN for CHUSJ preceded by '0'?
        else:
            patient = pho.get_patient_by_mrn(mrn)
        if patient is not None:
            pid = patient['id']
        else:
            logging.warning(f"Could not get PID using EP+MRN: {ep_mrn} nor by MRN: {mrn}.")
            # Retrieve PID using ramq?

    try:
        hpos = pho.parse_hpo(patient)
    except TypeError as e:
        logging.error(f"Could not use {ep_mrn} to retieve Phenotips patient: {patient}")
    else:
        for hpo in hpos:
            hpo_ids.append(hpo['id'])
            hpo_labels.append(hpo['label'])

    if len(hpo_ids) == 0:
        warn_msg = f"Could not find HPO terms for PID={pid} (EP+MRN={ep_mrn})"
        logging.warning(warn_msg)
        ids_str    = warn_msg
        labels_str = warn_msg
        logging.debug(f"Got HPO terms from Phenotips by Labeled EID {ep_mrn}\n")
        logging.debug(f"Phenotips ID for {ep_mrn} is {pid}")
        logging.debug(f"HPO labels_str is {labels_str}")
        logging.debug(f"HPO identifiers string is {ids_str}")
    else:
        ids_str = ';'.join(hpo_ids)
        labels_str = ';'.join(hpo_labels)

    return(pid, labels_str, ids_str)


def add_hpos_redcap(sample_name):
    """
    Lookup HPO identifiers from REDCap for `sample_name`.
    - sample_name : [str] Q1K sample name. Ex: 'Q1K_HSJ_10050_P'
    - Returns     : [str] Semi-column-spearated list of hpo identifiers
    """
    red  = REDCap()
    return red.get_hpo(sample_name)


def add_hpos_aoh():
    pass


def df_to_manifest(df):
    """
    From data in df, generate a manifest file for batch upload to Emedgene 
    using the UI. For specifications of the manifest, see:
    "https://help.emedgene.com/en/articles/7231644-csv-format-requirements"
    - `df`: A Pandas DataFrame
    - Returns: File 'emg_batch_manifest.csv' in current folder
    """
    df_manifest = pd.DataFrame({
        'Family Id': df['Family Id'],
        'Case Type': 'Whole Genome',
        'Files Names': df['filenames'],
        'Sample Type': 'FASTQ',
        'BioSample Name': df['sample_name'],
        'Visualization Files': '',
        'Storage Provider Id': 10126, # =prod. 10123=eval
        'Default Project': '',
        'Execute_now':  'False',
        'Relation': df['relation'],
        'Gender': df['gender'],
        'Phenotypes': 'Healthy',
        'Phenotypes Id': df['hpos'],
        #'Date Of Birth': pd.to_datetime(df['date_of_birth(YYYY-MM-DD)'], format='%d/%m/%Y'),
        'Date Of Birth': df['date_of_birth(YYYY-MM-DD)'],
        'Boost Genes': '',
        'Gene List Id': '',
        'Kit Id': '',
        'Selected Preset': 'Default',
        'Label Id': df['label'],
        'Clinical Notes': df['pid'],
        'Due Date': '',
        'Opt In': ''
    })
    # Convert Date of Birth to DateTime
    #
    try:
        df_manifest['Date Of Birth'] =  pd.to_datetime(df['date_of_birth(YYYY-MM-DD)'], format='%d/%m/%Y')
    except OverflowError as err:
        logging.warning(err)
    except:
        logging.warning(f"WARNING: Pandas could not convert Date of Birth column to DateTime format.")

    # With the "Files Names"="auto" option, BSSH users can automatically locate
    # FASTQ files based on the BioSample Name and Default Project provided.
    # Unfortunately, this would mean that cases woul bear the lab's CQGC_ID.
    #
    #df_manifest.loc[df_manifest['Relation'] == 'PROBAND', 'Phenotypes'] = df['phenotypes']
    df_manifest.loc[df_manifest['Relation'] == 'PROBAND', 'Phenotypes'] = ''
    df_manifest['Default Project'] = 'PRAGMatIQ_' + df_manifest['Label Id']

    df_manifest['Relation'].replace('PROBAND', 'proband', inplace=True)
    df_manifest['Relation'].replace('MTH', 'mother', inplace=True)
    df_manifest['Relation'].replace('FTH', 'father', inplace=True)
    df_manifest['Relation'].replace('BRO', 'sibling', inplace=True)
    df_manifest['Relation'].replace('SIB', 'sibling', inplace=True) # TODO: Verify 'SIB', or SIS?

    df_manifest['Gender'].replace('FEMALE', 'F', inplace=True)
    df_manifest['Gender'].replace('MALE', 'M', inplace=True)
    df_manifest['Gender'].replace('', 'U', inplace=True) 

    # Replace labels with corresponding IDs, which are platform-dependent
    # Use a correspondance table used to convert Labels to Label ID 
    # TODO: Use API to get list of codes instead of hard-coding the data
    #
    if args.site == 'prod':
        label2ID = {'CHUS': 12, 'CHUSJ': 13, 'CHUQ': 14, 'CUSM': 15}
    elif args.site == 'eval':
        label2ID = {'CHUS': 14, 'CHUSJ': 15, 'CHUQ': 16, 'CUSM': 17}
    else:
        logging.error(f"Option `--site|-s` ( '{args.site}') is not one of 'prod' or 'eval'")
    df_manifest['Label Id'] = df_manifest['Label Id'].apply(lambda x: label2ID[x])

    with open('emg_batch_manifest.csv', 'w') as fh:
        fh.write('[Data],,,,,,,,,,,,,,,,,,,,,\n')
        fh.write(df_manifest.to_csv(index=None, lineterminator='\n'))


def print_case_by_case(df):
    """
    Format and print df to STDOUT case by case, with HPO terms. 
    Easier reading, when creating cases manually using Emedgene's web UI.
    """
    pd.set_option('display.max_columns', 12)
    pd.set_option('display.max_colwidth', None)

    for case in df['Family Id'].unique():
        df_tmp = df[df['Family Id'] == case]
        pid    = df_tmp['pid'].tolist()[0]
        cohort = df_tmp['cohort_type'].tolist()[0]
        site   = df_tmp['label'].tolist()[0]
        print(f"============ {pid} | {case} | {site} | {cohort} ============\n")
        print(df_tmp[['pid', 'sample_name', 'biosample', 'relation', 'gender', 'date_of_birth(YYYY-MM-DD)', 'status']].to_string(index=False))
        hpo_terms = df_tmp[df_tmp['relation'] == 'PROBAND']['hpos']
        print(f"HPO Terms: {','.join(hpo_terms)}\n\n")


def main(args):
    """
    Read samples in --file and retrieve required information to build Cases:
    1. Add BaseSpace FASTQ file paths for each sample
    2. Get the corresponding HPO Identifiers and add HPO terms
    3. Use case PID instead of surname to sort and connect family members.
    4. Add storage provider and label ID's based on project
    Return a CSV file to be used as input for Emedgene's batch upload script.

    TODO: Check how QUADs are handled 
    TODO: Raise red flag when sibling or other family member is Affected 
    TODO: Add participants to cases
    TODO: Archive samples for this run
    """

    # Read samples in --file and retrieve required information to build Cases:
    #
    df_samples_list = pd.read_csv(args.file)
    workdir = os.path.dirname(os.path.abspath(args.file))
    os.chdir(workdir)
    print(f"\n# Log run {','.join(df_samples_list['flowcell'].unique())}\n")


    # 1. Add BaseSpace FASTQ file paths for each sample
    #
    logging.info(f"Add BaseSpace FASTQ file paths for each sample")
    df_samples_list['filenames'] = df_samples_list.apply(lambda row: add_fastqs(row.biosample), axis=1)
    logging.debug(f"Lookup FASTQ files for biosample {df_samples_list['biosample'][0]}:\n{add_fastqs(df_samples_list['biosample'][0])}")
    logging.debug(f"Filenames added as new column:\n{df_samples_list['filenames']}")


    # 2. Get the corresponding HPO Identifiers and add HPO terms
    #
    logging.info(f"Fetching HPO terms for project '{args.project}'")
    if args.project == 'prag' or args.project == 'eval':
        # HPO terms are stored in Phenotips for project PRAG. 
        # Also grab 'PID', which will populate 'Clinical Notes'
        #
        df_samples_list['hpos'] = df_samples_list.apply(lambda row: add_hpos_phenotips(row.ep_label, row.mrn)[2] if row.status == 'AFF' else '', axis=1)
        df_samples_list['pid']  = df_samples_list.apply(lambda row: add_hpos_phenotips(row.ep_label, row.mrn)[0] if row.status == 'AFF' else '', axis=1)
        logging.debug(f"Subset of DataFrame to show PIDs and HPO terms:\n{df_samples_list[['biosample', 'sample_name', 'status', 'pid', 'hpos']]}")
    elif args.project == 'q1k':
        # HPO terms are stored in REDCap for project Q1K.
        # add_hpos_redcap(sample_name) returns a semi-column-separated list of HPO terms.
        #
        df_samples_list['hpos'] = df_samples_list.apply(lambda row: add_hpos_redcap(row.sample_name) if row.status == 'AFF' else '', axis=1)
    elif args.project == 'aoh':
        # HPO terms are fixed.
        #
        pass
    else:
        logging.warn(f"Project '{args.project}' is not defined")


    # 3. Use case PID instead of surname to sort and connect family members.
    #
    logging.info(f"Sorting samples by families")


    # 4. Add storage provider and label ID's based on project
    #
    logging.info(f"Add storage provider and label ID's based for project {args.project}")
    
    # "Storage Provider ID" for BaseSpace and "Label ID" change according to 
    # the project domain (Organization) we're logged into. Dict to set the IDs
    # based on `args.project`
    # Maybe this should be in the configuration file?
    #
    projects_ids = {
        'eval': {'storage_id': '10123', 'label_ids': {'CHUS': '14', 'CHUSJ': '15', 'CHUQ': '16', 'CUSM': '17'}},
        'prag': {'storage_id': '10126', 'label_ids': {'CHUS': '12', 'CHUSJ': '13', 'CHUQ': '14', 'CUSM': '15'}},
        'q1k' : {'storage_id': '10219', 'label_ids': {'CHUSJ': '1'}},
        'aoh' : {'storage_id': '10220', 'label_ids': {'CHUSJ': '1'}}
    }

    # Return a CSV file to be used as input for Emedgene's batch upload script
    #
    logging.info("Wrote manifest file `emg_batch_manifest.csv` for batch upload to Emedgene.")
    logging.info(f"Done.\n")


def tests():
    return(1)

if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)
    #tests()

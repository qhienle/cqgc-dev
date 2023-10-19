#!/usr/bin/env python3
"""
Get Case information from Nanuq with Phenotips (PID) for a given _Run_.

USAGE: get_nanuq_files.py A00516_420
       get_nanuq_files.py --help

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
from lib.gapp import Phenotips
from lib.gapp import BSSH

nq = Nanuq()
#pho = Phenotips()
bssh = BSSH()

__version__ = "0.2"


def parse_args():
    parser = argparse.ArgumentParser(description="Get Case information from Nanuq for a given Run.")
    parser.add_argument('run', help="FC_SHORT Run ID, ex: 'A00516_339'")
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


def format_mrn_eid(ep, mrn):
    """
    Format Medical Record Number (MRN) identifiers for Phenotips
    - ep (Etablissement Public): Site label [str]
    - mrn: Medical Record Number [str]
    - Returns: Phenotips identifier [str]

    MRN is used as identifier to look up a record in Phenotips (by API method 
    labeled external ID, "labeled-eid"). In Phenotips, MRN is prepended with 
    the Site's EP initials (_e.g._ "CHUS1626861"). There does not seem to be a
    standard format for MRN identifiers. Here are some of the format detected, 
    by order of frequency of occurrences:

    PRAGMatIQ | Phenotips     | Nanuq       | Notes
    ----------|---------------|-------------|------------------------------------
    3421069   | CHUSJ3421069  | 03421069    | For CHUSJ, Nanuq adds a leading "0"
    X3627954  | CHUSJX3627954 | X3627954    | Not numerical, starts with 'X'
    1628699   | CHUS1628699   | 1628699     | 7 digits, no leading '0'
    347990    | CHUS347990    | 347990      | 6 and no leading '0' added by Nanuq
    1633799   | 1633799       | 1633799     | Not prefixed with EP initials
    1644460   | CHUQ1644460   | CHUL1644460 | CHUQ is stored as CHUL in Nanuq
    """
    if ep == 'CHUSJ':
        mrn = mrn.lstrip('0')
    elif ep == 'CHUS':
        pass
    elif ep == 'CHUQ':
        return(mrn.replace('L', 'Q'))
    return(ep + mrn)


def add_hpos(ep, mrn):
    """
    Lookup Phenotips ID (PID) and HPO identifiers, using `ep_mrn`
    - ep_mrn : [str] EP+MRN identifier. Ex: CHUSJ123456
    - Returns:
    """
    pho = Phenotips()
    pid = ''
    hpo_ids    = []
    hpo_labels = []

    ep_mrn   = format_mrn_eid(ep, mrn)
    patient  = pho.get_patient_by_mrn(ep_mrn)
    warn_msg = f"Could not get PID using EP+MRN: {ep_mrn}"

    if patient is not None:
        pid = patient['id']
        hpos = pho.parse_hpo(patient)
        for hpo in hpos:
            hpo_ids.append(hpo['id'])
            hpo_labels.append(hpo['label'])
    else:
        logging.warning(warn_msg)

    if len(hpo_ids) == 0:
        ids_str    = warn_msg
        labels_str = warn_msg
    else:
        ids_str = ','.join(hpo_ids)
        labels_str = ','.join(hpo_labels)

    # sample_infos.append(pid)
    # sample_infos.append(labels_str)
    # sample_infos.append(ids_str)
    logging.debug(f"Got HPO terms from Phenotips by Labeled EID {ep_mrn}\n")
    logging.debug(f"Phenotips ID for {ep_mrn} is {pid}")
    logging.debug(f"HPO labels_str is {labels_str}")
    logging.debug(f"HPO identifiers string is {ids_str}")
    
    return(pid, labels_str, ids_str)


def df_to_manifest(df):
    """
    Generate a manifest file for batch upload to Emedgene from data in df.
    See "Case_creation-script_v2.docx" for manifest specifications.
    - `df`: A Pandas DataFrame
    - Returns: File 'emg_batch_manifest.csv' in current folder
    """
    df_manifest = pd.DataFrame(
        {
            'case_group_number': df['case_group_number'],
            'case_type': 'Whole Genome',
            'filenames': df['filenames'],
            'bam_file': '', 
            'execute_now':  'False',
            'sample_name': df['sample_name'],
            'relation': df['relation'],
            'gender': df['gender'],
            'phenotypes': df['phenotypes'],
            'hpos': df['hpos'],
            'boost_genes': '',
            'gene_list_id': '',
            'kit_id': '',
            'selected_preset': '',
            'due_date(YYYY-MM-DD)': '',
            'label': df['label'],
            'bigwig': '',
            'clinical_notes': df['pid'],
            'Default Project': '',
            'date_of_birth(YYYY-MM-DD)': df['date_of_birth(YYYY-MM-DD)'],
        }
    )
    df_manifest['Default Project'] = 'PRAGMatIQ_' + df_manifest['label']

    with open('emg_batch_manifest.csv', 'w') as fh:
        fh.write(df_manifest.to_csv(index=None, lineterminator='\n'))
    
    
def print_case_by_case(df):
    """
    Format and print df to STDOUT case by case, with HPO terms. 
    Easier reading, when creating cases manually using Emedgene's web UI.
    """
    pd.set_option('display.max_columns', 12)
    pd.set_option('display.max_colwidth', None)

    for case in df['case_group_number'].unique():
        df_tmp = df[df['case_group_number'] == case]
        pid    = df_tmp['pid'].tolist()[0]
        cohort = df_tmp['cohort_type'].tolist()[0]
        site   = df_tmp['label'].tolist()[0]
        print(f"============ Case INFO: {case} | {site} | {cohort} | {pid} ============\n")
        print(df_tmp[['pid', 'sample_name', 'biosample', 'relation', 'gender', 'date_of_birth(YYYY-MM-DD)', 'status']].to_string(index=False))
        hpo_terms = df_tmp[df_tmp['relation'] == 'PROBAND']['hpos']
        print(f"\nHPO Terms:\n{','.join(hpo_terms)}\n\n")


def list_samples_to_archive(df):
    """
    Create a file listing samples to archive. This list can later be used to
    archive samples and to collect metrics with scripts `archive_PRAGMatIQ.sh`
    and `emg_collect_samples_metrics.py`, respectively.
    - `df`: A Pandas DataFrame
    - Returns: list of samples [str] and file 'samples_list.txt' for archiving
    """
    filename = 'samples_list.csv'
    df1 = df[['sample_name', 'biosample', 'label', 'fc_date']] # TODO: Add flowcell
    df1 = df1.rename(columns={'sample_name': 'Sample', 'biosample': 'CQGC_ID', 'label': 'Site', 'fc_date': 'Date'})
    df1.to_csv(filename, index=False)
    logging.info(f"Created file {filename}")
    return(f"{' '.join(df['Sample'])}")


def main(args):
    """
    Retrieve necessary information from Nanuq for creating cases in Emedgene.
    Return a CSV file to be used as input for Emedgene's batch upload script.

    1. Download list of samples from Nanuq (API) for a given Run ID.
    2. For each sample in the list of SampleNames: 
        2.1 Get the JSON from Nanuq to extract infos to create cases on EMG;
        2.2 Get the Phenotips ID (PID) and the corresponding HPO Identifiers;
        2.3 Connect to BaseSpace and re-construct the path to the FASTQ files;
    3. Combine individual data into a Pandas data frame
        3.1 Sort, group and print each trio to STDOUT for case creation.
        3.2 Use case PID instead of surname to connect family members.
    4. Convert DataFrame into a CSV file (manifest) for EMG batch upload.
    5. TODO: Batch upload to Emedgene using their script
    6. TODO: Add participants to cases
    7. TODO: Archive samples for this run
    """

    print(f"# Logging run {args.run}")

    workdir = f"{os.getcwd()}{os.sep}{args.run}"
    try:
        os.mkdir(workdir)
        logging.info(f"Created work directory '{workdir}'")
    except:
        logging.warning(f"Could not create {workdir}")
    os.chdir(workdir)

    # 1. Get a list of samples on this run to construct the cases.
    # TODO: Add experiment name as an alternative identifier for Nanuq API?
    #
    samplenames = nq.get_samplenames(args.run)
    if not samplenames.text.startswith("##20"):
        sys.exit(logging.error(f"Unexpected content for SampleNames. Please verify Nanuq's reponse:\n{samplenames.text}"))
    else:
        logging.info("Retrieved samples conversion table from Nanuq")

    fc_date = re.match(r'##(\d{4}-\d{2}-\d{2})', samplenames.text).group(1)
    logging.debug(f"Date of run from Nanuq's SampleNames file: {fc_date}")

    
    # 2. Build cases: Get Nanuq JSON for each CQGC ID found in SampleNames 
    # (returned as a string by requests.text) and parse sample infos. 
    # SampleNames lines are tab-delimitted. Comment lines begin with "#".
    # Results are stored in `cases` and printed to STDOUT at the end.
    # 
    cases = []
    for line in samplenames.text.splitlines():
        if not line.startswith('#'):
            cqgc, sample = line.split("\t")
            
            # 2.1 Get information for sample from Nanuq
            #
            try:
                data = json.loads(nq.get_sample(cqgc))
            except Exception as e:
                logging.warning(f"JSONDecodeError {e} could not decode sample {cqgc} ({sample})")
                continue

            logging.info(f"Got information for biosample {cqgc} a.k.a. {sample}")
            if len(data) != 1:
                logging.debug(f"Number of samples retrieved from Nanuq is not 1.\n{data}")
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

            # 2.2 Add Phenotips ID (`pid`) and patients' HPO identifiers for
            # the proband. Lookup this information in Phenotips, using EP+MRN
            # Ex: CHUSJ123456
            #
            if data[0]["patient"]["familyMember"] == 'PROBAND':
                pid, labels_str, ids_str = add_hpos(data[0]["patient"]["ep"], data[0]["patient"]["mrn"])
                logging.info(f"Got HPO terms from Phenotips for PID {pid}")
            else:
                pid, labels_str, ids_str = ('', '', '')
                logging.debug(f'Not retrieving PID for {cqgc} ({data[0]["patient"]["familyMember"]})')
            sample_infos.append(pid)
            sample_infos.append(labels_str)
            sample_infos.append(ids_str)
            logging.debug(f"PID: {pid}; HPO ID: {ids_str}; Labels: {labels_str}\n")

            # 2.3 Add paths to fastq on BaseSpace
            #
            fastqs = bssh.get_sequenced_files(data[0]["labAliquotId"])
            sample_infos.append(';'.join(fastqs))

            cases.append(sample_infos)
    
    # 3. Load cases (list of list) in a DataFrame, sort and group members
    # Translate column names to match EMG's manifest specifications.
    # pid => case_group_number, hpo_labels => phenotypes, hpo_ids => hpos
    # Group by family and sort by relation.
    # Add case_group_number (PID) to all family members based on familyID.
    #
    df = pd.DataFrame(cases)
    df.columns = ['sample_name', 'biosample', 'relation', 'gender', 'label', 
                  'mrn', 'cohort_type', 'date_of_birth(YYYY-MM-DD)', 'status',
                  'case_group_number', 'pid', 'phenotypes', 'hpos', 'filenames']
    df['fc_date'] = fc_date
    logging.info(f"Add column for flowcell date {fc_date}")
    df = df.sort_values(by=['case_group_number', 'relation'], ascending=[True, False])
    #logging.info("Sorted families. Setting PID as case_group_number")
    
    # Print to STDOUT case by case, with HPO terms. Easier reading, when 
    # creating cases manually using Emedgene's web UI
    #
    logging.info(f"\nCases for {args.run}:\n")
    df1 = df.drop(['phenotypes', 'filenames'], axis=1)
    print_case_by_case(df1)

    # 4. Output manifest for batch upload, see "Case_creation-script_v2.docx"
    #
    df_to_manifest(df)

    # 5. Batch upload to Emedgene using their script
    #
    logging.info("Please run the command below, replacing '-u USER' and '-p PASS' with Emedgene credentials:")
    print('python /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/create_batch_cases_v2.py -i emg_batch_manifest.csv -s 10123 -hu stejustine.emedgene.com -u cqgc.bioinfo.hsj@ssss.gouv.qc.ca -p PASS -b\n')
    # subprocess.run(['python', '/staging2/soft/CQGC-utils/Analysis.pipeline_illumina/create_batch_cases_v2.py', 
    #                 '-i', 'emg_batch_manifest.csv', 
    #                 '-s', '10123', 
    #                 '-hu', 'stejustine.emedgene.com', 
    #                 '-u', 'cqgc.bioinfo.hsj@ssss.gouv.qc.ca', 
    #                 '-p', 'PASS', 
    #                 '-b'])

    # TODO: 6. Add participants to cases
    #
    
    # TODO: 7. Archive samples from cases finalized on Emedgene
    #
    logging.info(f"List of samples to archive:\n{list_samples_to_archive(df1)}")


def tests():
    print("Running in test mode")

    # 1. Get a list of samples on this run to construct the cases.
    # TODO: Add experiment name as an alternative identifier for Nanuq API?
    #
    samplenames = nq.get_samplenames(args.run)
    return(samplenames)

if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)

    main(args)
    #tests()


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
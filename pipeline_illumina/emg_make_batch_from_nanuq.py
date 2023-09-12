#!/usr/bin/env python3
"""
Get Case information from Nanuq with Phenotips (PID) for a given _Run_.

USAGE: get_nanuq_files.py --run A00516_420
       get_nanuq_files.py --help

Nanuq username and password have be saved in a file named '~/.nanuq', like so:
`echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`
Replace USERNAME and PASSWORD with actual values.

Tokens to connect with Phenotips and BaseSpace (BSSH) are expected to be found 
in ~/.illumina/gapp_conf.json (available at https://github.com/CQGC-Ste-Justine/PrivateDoc/)
"""

import os, sys
import time
import argparse
import subprocess
import json
import pandas as pd

# Set source path to CQGC-utils so that we can use relative imports
#
src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(src_path)
from lib.nanuq import Nanuq
from lib.gapp import Phenotips
from lib.gapp import BSSH

__version__ = "0.2"


def parse_args():
    parser = argparse.ArgumentParser(description="Get Case information from Nanuq for a given Run.")
    parser.add_argument('-r', '--run', required=True, help="FC_SHORT Run ID, ex: 'A00516_339'")
    return(parser.parse_args())


def now():
    """
    Return Date-Time string for logging. Ex.: print(f"{now()} Hello World!").
    """
    return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))


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
    by order of frequency of occurences:

    PRAGMatIQ | Phenotips     | Nanuq    | Notes
    ----------|---------------|----------|------------------------------------
    3421069   | CHUSJ3421069  | 03421069 | For CHUSJ, Nanuq adds a leading "0"
    X3627954  | CHUSJX3627954 | X3627954 | Not numerical, starts with 'X'
    1628699   | CHUS1628699   | 1628699  | 7 digits, no leading '0'
    347990    | CHUS347990    | 347990   | 6 and no leading '0' added by Nanuq
    1633799   | 1633799       | 1633799  | Not prefixed with EP initials
    """
    if ep == 'CHUSJ':
        mrn = mrn.lstrip('0')
    elif ep == 'CHUS':
        ep = 'CHUQ'
    return(ep + mrn)


def main(args):
    """
    Retrieve necessary information from Nanuq for creating cases in Emedgene.
    Return a CSV file to be used as input for Emedgene's batch upload script.

    1. Download list of samples from Nanuq (API) for a given Run ID.
    2. For each sample in the list of SampleNames: 
        2.1 Get the JSON from Nanuq to extract infos to create cases on EMG;
        2.2 Get the Phenotips ID (PID) and the corresponding HPO Identifiers;
        2.3 Associate patient surname to PID, for later use
        2.4 Connect to BaseSpace and re-construct the path to the FASTQ files;
    3. Combine individual data into a Pandas dataframe
        3.1 Sort, group and print each trio to STDOUT for case creation.
        3.2 Use case PID instead of surname to connect family members.
    4. Convert DataFrame into a CSV file (manifest) for EMG batch upload.
    5. TODO: Batch upload to Emedgene using their script
    6. TODO: Add participants to cases
    7. TODO: Archive samples for this run
    """

    bssh = BSSH()      # Handler to work with BSSH
    nq   = Nanuq()     # Interact with Nanuq REST API
    pho  = Phenotips() # Interact with Phenotips REST API

    # PID is used to group family members, instead of the family name
    # (nominative info). Build a lookup table to assign PID to members with
    # the same family name.
    #
    familyId2pid = {}   # Lookup table {'surname': 'pid', 'surname': 'pid',...}

    # 1. Get a list of samples on this run to construct the cases.
    # TODO: Add experiment name as an alternative identifier for Nanuq API?
    #
    samplenames = nq.get_samplenames(args.run)
    if not samplenames.text.startswith("##20"):
        sys.exit(f"{now()} ERROR: Unexpected content for SampleNames. Please verify Nanuq's reponse:\n{samplenames.text}")
    else:
        print(f"{now()} Retrieved samples conversion table from Nanuq")
    
    # 2. Build cases: Get Nanuq JSON for each CQGC ID found in SampleNames 
    # (returned as a string by requests.text) and parse sample infos. 
    # SampleNames lines are tab-delimitted. Comment lines begin with "#".
    # Results are stored in `cases` and printed to STDOUT at the end.
    # 
    cases = []
    for line in samplenames.text.splitlines():
        if not line.startswith('#'):
            cqgc, sample = line.split("\t")
            
            # 2.1 Get information for sample frm Nanuq
            #
            data = json.loads(nq.get_sample(cqgc))
            print(f"{now()} Got information for biosample {cqgc} a.k.a. {sample}")
            if len(data) != 1:
                print(f"{now()} WARNING: Number of samples retrieved from Nanuq is not 1.\n{data}")
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

            # 2.2 Add Phenotips ID (`pid`) and patients' HPO identifiers
            # Lookup this information in Phenotips, using the EP+MRN
            # Ex: CHUSJ123456
            #
            ep_mrn = format_mrn_eid(data[0]["patient"]["ep"], data[0]["patient"]["mrn"])
            print(f"{now()} Getting HPO terms from Phenotips by Labeled EID {ep_mrn}\n")

            patient   = pho.get_patient_by_mrn(ep_mrn)
            hpo_ids   = []
            hpo_labels= []
            if patient is not None:
                pid = patient['id']
                hpos = pho.parse_hpo(patient)
                for hpo in hpos:
                    hpo_ids.append(hpo['id'])
                    hpo_labels.append(hpo['label'])
            else:
                pid = ''

            if len(hpo_ids) == 0:
                ids_str = ''
                labels_str= ''
            else:
                ids_str = ','.join(hpo_ids)
                labels_str = ','.join(hpo_labels)

            sample_infos.append(pid)
            sample_infos.append(labels_str)
            sample_infos.append(ids_str)

            # 2.3 Add family name and PID to the lookup table
            #
            familyId = data[0]['patient']['familyId']
            if familyId not in familyId2pid and pid.startswith('P'):
                familyId2pid[familyId] = pid

            # 2.4 Add paths to fastq on BaseSpace
            #
            fastqs = bssh.get_sequenced_files(data[0]["labAliquotId"])
            sample_infos.append(';'.join(fastqs))

            cases.append(sample_infos)

    # 3. Load cases (list of list) in a DataFrame, sort and group members
    # Translate column names to match EMG's manifest specifications.
    # Group by family and sort by relation.
    # Add case_group_number (PID) to all family members based on familyID.
    #
    df  = cases_to_df(cases)
    df = pd.DataFrame(cases)
    df.columns = ['sample_name', 'biosample', 'relation', 'gender', 'label', 
                  'mrn', 'cohort_type', 'date_of_birth(YYYY-MM-DD)', 'status',
                  'family', 'case_group_number', 'phenotypes', 'hpos', 'filenames']
    df = df.sort_values(by=['family', 'relation'], ascending=[True, False])
    print(f"{now()} Sorted families. Set PID as case_group_number based on look up table familyId2pid:\n{familyId2pid}")
    for index, row in df.iterrows():
        if row['case_group_number'] == '':
            try:
                row['case_group_number'] = familyId2pid[row['family']]
            except KeyError as err:
                print(f"{now()} ***WARNING!*** Could not set PID as family identifier. KeyError: {err}")
    # try:
    #     df['case_group_number'] = df.apply(lambda x: familyId2pid[x['family']], axis=1)
    # except KeyError as err:
    #     print(f"{now()} ***WARNING!*** Could not set PID as family identifier. KeyError: {err}")
    # else:
    #     print(f"{now()} Sorted families and assigned case_group_number based on PID\n{familyId2pid}")

    pd.set_option('display.max_columns', 12)
    pd.set_option('display.max_colwidth', None)

    df1 = df.drop(['phenotypes', 'filenames'], axis=1)
    print(f"\n{now()} Cases for {args.run}:\n")
    # print(df.drop(['phenotypes', 'filenames'], axis=1))
    

    
    print(f"{now()} List of samples to archive after cases are finalized on Emedgene: {df1['sample_name']}")

    # 4. Output manifest for batch upload, see "Case_creation-script_v2.docx"
    #
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
            'clinical_notes': df['case_group_number'],
            'Default Project': '',
            'date_of_birth(YYYY-MM-DD)': df['date_of_birth(YYYY-MM-DD)'],
        }
    )
    df_manifest['Default Project'] = 'PRAGMatIQ_' + df_manifest['label']

    with open('emg_batch_manifest.csv', 'w') as fh:
        fh.write(df_manifest.to_csv(index=None, lineterminator='\n'))

    # 5. Batch upload to Emedgene using their script
    #
    print(f"{now()} Please run the command below, replacing '-u USER' and '-p PASS' with Emedgene credentials:")
    print('\npython /staging2/soft/CQGC-utils/Analysis.pipeline_illumina/create_batch_cases_v2.py\\') 
    print('\t-i emg_batch_manifest.csv \\')
    print('\t-s 10123 \\')
    print('\t-hu stejustine.emedgene.com \\')
    print('\t-u cqgc.bioinfo.hsj@ssss.gouv.qc.ca \\')
    print('\t-p PASS \\')
    print('\t-b\n')
    # subprocess.run(['python', 
    #                 '/staging2/soft/CQGC-utils/Analysis.pipeline_illumina/create_batch_cases_v2.py', 
    #                 '-i', 'emg_batch_manifest.csv', 
    #                 '-s', '10123', 
    #                 '-hu', 'stejustine.emedgene.com', 
    #                 '-u', 'cqgc.bioinfo.hsj@ssss.gouv.qc.ca', 
    #                 '-p', 'PASS', 
    #                 '-b'])

    # TODO: 6. Add participants to cases
    
    # TODO: 7. Archive samples from cases finalized on Emedgene
    #
    print(f"{now()} Please run the command below on narval to archive samples from finalized csaes on this run:")
    print(f"{df1['sample_name']}")
    
    
def tests(args):
    print("Running in test mode")
    nq = Nanuq()
    print(nq.get_sample(21571))
    print("\nDone.\n")


if __name__ == '__main__':
    args = parse_args()
    main(args)
    #tests(args)

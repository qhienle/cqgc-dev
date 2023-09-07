import base64
from datetime import datetime
import os
from logging.handlers import RotatingFileHandler

import requests
import argparse
import csv
import json
import pathlib
import logging
import sys
from itertools import groupby

URL = ''
RELATIONS = ['proband', 'father', 'mother', 'sibling']
is_basespace_storage_provider = False
phenotype_extractor = None
phenotype_extractor_type = None
basespace_biosample_name_to_sequenced_files_map = {}

VERSION = "2.6"

DEBUG_VERSION = False
GENDER_DICT = {'m': 'Male', 'f': 'Female', 'u': 'Unknown'}

ROUTE_AUTH_TOKEN = 'https://{}/api/auth/api_login/'
ROUTE_USERS = 'https://{}/api/organization/users/'
ROUTE_GET_STORAGE_RESOURCES = 'https://{}/api/storage/{}/fs/list'
ROUTE_CREATE_CASE = 'https://{}/api/cases/v2/cases'
ROUTE_GET_PHENOTYPE = 'https://{}/api/phenotype/?query={}'

CASE_GROUP_NUMBER = 'case_group_number'
EXECUTE_NOW = 'execute_now'
CASE_TYPE = 'case_type'
FILES_NAMES = 'files_names'
BAM_FILE = 'bam_file'
SAMPLE_NAME = 'sample_name'
DEFAULT_PROJECT = 'Default Project'
RELATION = 'relation'
GENDER = 'gender'
PHENOTYPES = 'phenotypes'
HPOS = 'hpos'
BOOST_GENES = 'boost_genes'
GENE_LIST_ID = 'gene_list_id'
KIT_ID = 'kit_id'
SELECTED_PRESET = 'selected_preset'
DUE_DATE = 'due_date(YYYY-MM-DD)'
LABELS = 'labels'
BIGWIG = 'bigwig'
CLINICAL_NOTES = 'clinical_notes'
DATE_OF_BIRTH = 'date_of_birth(YYYY-MM-DD)'


# (*) - mandatory field
# header fields - 20 fields
##########################################
# case_group_number (*) - starts from 1 and increment by 1 for the next case group. for pedigree - set the same ordinal number
# execute_now - True | False - default is True, False will set case into pending-sequencing !
# case_type (*) - Exome | Whole Genome | Custom Panel
# files_names - list of sample file paths, semicolon (;) separated
# bam_file - path to a bam file (for vcf sample type)
# sample_name (*) - sample name
# Default Project - project name (basespace compatibility)
# relation (*) - proband | father | mother | sibling
# gender (*) - M | F | U
# phenotypes - unaffected | phenotypes name list, semicolon (;) separated. (you may use only hpos)
# hpos - hpo id list, semicolon (;) separated. (you may use only hpos)
# boost_genes - True | False - default False
# gene_list_id - id of gene list - need to get the id from customer support - default is all genes
# kit_id - id of a kit - need to get the id from customer support
# selected_preset - id of the selected preset - need to get the id from customer support
# due_date(YYYY-MM-DD) - case due date format: YYYY-MM-DD (2022-12-30)
# labels - list of labels IDs, semicolon (;) separated. (for proband only)
# bigwig - path to a bigwig file
# clinical_notes - clinical notes
# date_of_birth(YYYY-MM-DD) - patient date of birth. format: YYYY-MM-DD (2022-12-30)

# all paths are related to the storage provider that is given in the command line arguments to this script
# command line example:
# -hu {host_name}.emedgene.com -u {username} -p {password} -s {storage_provider_id} -i {metadata_file_name}.csv|json

# if using 'Phenotips' to get the phenotypes you should add to the command line the following arguments:
# -t1 {phenotypes_server_type} -d1 {phenotypes_domain} -u1 {phenotypes_username} -p1 {phenotypes_password}
# -s1 {phenotypes_secret}
# and set the hpo column with the eid


def parse_args():
    parser = argparse.ArgumentParser(description='This script is used for creating a batch of cases '
                                                 'using csv/json that contains patients data. Script supports '
                                                 'Basespace and all other storage provider that are supported by '
                                                 'the Emedgene Platform')
    parser.add_argument('-i', '--info-file', dest='info_file_path', required=True,
                        help='Path to csv/json file with case information.')
    parser.add_argument('-s', '--storage-id', dest='storage_id', required=True,
                        help='ID of the storage provider that hold the files.')
    parser.add_argument('-hu', '--host-url', dest='host_url', required=False, default=URL,
                        help='Emedgene host.')
    parser.add_argument('-u', '--username', dest='username', required=True,
                        help='Username in platform.')
    parser.add_argument('-p', '--password', dest='password', required=True,
                        help='User password matches the username.')
    parser.add_argument('-b', '--is_basespace', dest='is_basespace', action='store_true', required=False, default=False,
                        help='If flag exists: Use basespace convention, else use standard convention')
    parser.add_argument('-v', '--version', dest='version', action='store_true', help='version')

    # phenotypes added separately
    parser.add_argument('-t1', '--phenotype_extractor_type', dest='phenotype_extractor_type', required=False,
                        default=None)
    parser.add_argument('-d1', '--phenotypes-domain', dest='phenotypes_domain', required=False, default='')
    parser.add_argument('-u1', '--phenotypes-username', dest='phenotypes_username', required=False, default='')
    parser.add_argument('-p1', '--phenotypes-password', dest='phenotypes_password', required=False, default='')
    parser.add_argument('-s1', '--phenotypes-secret', dest='phenotypes_secret', required=False, default='')

    return parser.parse_args()


def setup_logger(log_folder, log_prefix):
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    for handler in root.handlers:
        root.removeHandler(handler)
    # create logfile handler
    if not os.path.exists(log_folder):
        os.mkdir(log_folder)
    log_name = f'{log_folder}/{log_prefix}-{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.log'
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logfile_handler = RotatingFileHandler(log_name, maxBytes=50000000, backupCount=5)
    logfile_handler.setLevel(logging.INFO)
    logfile_handler.setFormatter(formatter)
    root.addHandler(logfile_handler)

    # create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)


class PhenotypeExecutor:
    # this class is only suitable for 'Phenotips' at the moment
    ROUTE_GET_PATIENT_DATA_BY_LABEL_EID = 'https://{}/rest/patients/labeled-eid/{}/{}'

    def __init__(self, extractor_type, domain, username, password, secret):
        self.extractor_type = extractor_type
        self.domain = domain
        self.user_name = username
        self.password = password
        self.secret = secret
        self.label = "DNA number"
        self.eid = ""
        self.headers = self.get_user_authentication_headers()

    def get_user_authentication_headers(self):
        # prepare the encoded username & password
        username_password = f'{self.user_name}:{self.password}'
        username_password_bytes = username_password.encode("ascii")
        username_password_base64_bytes_encoded = base64.b64encode(username_password_bytes)
        username_password_base64_string_encoded = username_password_base64_bytes_encoded.decode("utf-8")

        headers = {
            "accept": "application/json",
            "authorization": f"Basic {username_password_base64_string_encoded}",
            "X-Gene42-Secret": self.secret
        }
        return headers

    def get_hpo_str(self):
        url = self.ROUTE_GET_PATIENT_DATA_BY_LABEL_EID.format(self.domain, self.label, self.eid)
        response = requests.get(url=url, headers=self.headers)
        if response.status_code != 200:
            logging.error(f"get_hpo_str failed response code: {response.status_code}, message: {response.text}")
            response.raise_for_status()

        return ";".join([elem.get("id") for elem in response.json().get("features")])

    def set_phenotype_eid(self, eid):
        self.eid = eid


def get_authentication_token(username, password, host_url):
    logging.info('Getting the Authorization bearer token')
    url = ROUTE_AUTH_TOKEN.format(host_url)
    payload = {"username": username, "password": password}
    response = requests.post(url, json=payload)
    response.raise_for_status()
    auth_header = response.json().get("Authorization")
    return auth_header


def get_remote_object_list(url, auth_token, path):
    has_more = True
    from_value = 0
    size_value = 100
    resources_list = []
    while has_more:
        params = {"from": from_value, "size": size_value, "path": path}
        response = requests.get(url, params=params, headers={'Authorization': auth_token})
        if response.status_code != 200:
            logging.error(f"get_remote_file_list failed response code: {response.status_code} message: {response.text}")
            response.raise_for_status()
        else:
            logging.info(f"get_remote_file_list passed - path: {path}")
            resources_list.extend(response.json().get('resources'))
            has_more = response.json().get('has_more')
            if has_more:
                from_value += size_value

    return resources_list


def recursive_construct_basespace_path_with_ids(current: list, split_path: list, storage_url: str, auth_token: str):
    if not split_path:
        # nothing lefty to check. return the current path
        return current
    base = split_path[0]
    end = split_path[1:]
    response_file_list = get_remote_object_list(storage_url, auth_token, "/".join(current))
    for file in response_file_list:
        if file.get("alias").lower() == base.lower():
            current.append(file.get("name"))
            final_path = recursive_construct_basespace_path_with_ids(current, end, storage_url, auth_token)
            if final_path:
                # the current path found a match on BSSH.
                return final_path
    # this path didn't find a match on BSSH, so we need to pop it from the current
    # this is because the list is mutable and is used in all the different recursive calls

    if len(current) == 0:
        raise Exception(f"path not found on basespace. path: {'/'.join(split_path)}")

    current.pop()
    return


def convert_basespace_path_to_path_with_ids(host_url, storage_id, auth_token, human_readable_path):
    storage_url = ROUTE_GET_STORAGE_RESOURCES.format(host_url, storage_id)
    split_path = [step for step in human_readable_path.split('/') if step.strip() != ""]
    id_path_list = recursive_construct_basespace_path_with_ids([], split_path, storage_url, auth_token)
    logging.info(f'{human_readable_path} -> {"/" + "/".join(id_path_list)}')
    return "/" + "/".join(id_path_list)


def extract_phenotypes(item_list, host_url, auth_token):
    phenotype_hpo_list = []
    for item in item_list:
        term = item.upper().strip()
        phenotype_url = ROUTE_GET_PHENOTYPE.format(host_url, term)
        response = requests.get(phenotype_url, headers={'Authorization': auth_token})
        if response.status_code != 200:
            logging.error(f"extract_phenotypes failed response code: {response.status_code} message: {response.text}")
            response.raise_for_status()
        else:
            res = response.json()
            count = res.get('total')
            if count == 0:
                error = f"Extract_phenotypes for HPO/Phenotype: '{item}' returned count: {count} phenotypes-hpo"
                raise ValueError(error)
            elif count > 1:
                for disease in res.get('hits'):
                    if disease.get('name').upper() == term:
                        if 'match' in disease.keys():
                            del disease['match']
                        phenotype_hpo_list.append(disease)
                        break
            else:  # only 1 item
                disease = res.get('hits')[0]
                if 'match' in disease.keys():
                    del disease['match']
                phenotype_hpo_list.append(disease)

    return phenotype_hpo_list


def get_phenotype_objects(family_member, host_url, auth_token):
    if phenotype_extractor:
        eid = family_member.get(HPOS)
        phenotype_extractor.set_phenotype_eid(eid)
        hpos = phenotype_extractor.get_hpo_str()
    else:
        phenotype_names = family_member.get(PHENOTYPES)
        hpos = family_member.get(HPOS)
    # get phenotypes from phenotypes api
    if hpos:
        return extract_phenotypes(hpos.split(';'), host_url, auth_token)
    # no hpos
    elif phenotype_names and phenotype_names.lower() != 'unaffected':
        return extract_phenotypes(phenotype_names.split(';'), host_url, auth_token)
    return []


def get_ethnicity():
    return {'maternal': [], 'paternal': []}


def get_base_json():
    return {
        "analysis_type": None,  # (None, 'carrier', 'rare_disease', 'population_screening', 'from_kit', 'virtual_panel')
        "boostGenes": False,
        "consanguinity": False,
        # "disease_penetrance": 0,
        "disease_severity": "",
        # "diseases": [],
        # "due_date": "",
        "extra_data": "",
        "filename": "samples/",
        "gene_list": {
            "type": "all",
            "id": 1,
            "visible": False
        },
        # "gene_list_visible": False,
        # "incidental_findings": "Yes",
        "indication_for_testing": "string",
        "inheritance_modes": [],
        "labels": [],
        "notes": "",
        "patients": {
            "other": []
        },
        "sample_type": "",
        "samples": [],
        "selected_preset_set": "Default",
        "sequence_info": {},
        "type": ""
    }


def set_sample_type(base_json, files_names, case_group_number):
    sample_file_name = files_names.split(";")[0]
    if not sample_file_name:
        logging.info(f"Case group number: [{case_group_number}] - no samples files are set")
    elif sample_file_name.endswith('fastq.gz') or sample_file_name.endswith('fq.gz'):
        base_json['sample_type'] = "fastq"
        logging.info(f"Case group number: [{case_group_number}] - sample type: fastq")
    elif sample_file_name.endswith('vcf.gz') or sample_file_name.endswith('vcf'):
        base_json['sample_type'] = "vcf"
        logging.info(f"Case group number: [{case_group_number}] - sample type: vcf")


def get_project_id(family_member, storage_url, auth_token, project_list):
    project_name = family_member.get(DEFAULT_PROJECT, '')
    sample_name = family_member.get(SAMPLE_NAME, '').strip()

    # loop over all projects in order to find the project with the given project id
    for project in project_list:
        if project_name:  # if given the project name return the project id
            if project_name.lower() == project.get('alias').lower():
                return project.get('name')
        else:
            # search for the project, where the given sample is, and return that project id
            sample_list = get_sample_list(storage_url, auth_token, project.get('name'))
            for sample in sample_list:
                if sample.get('alias').lower() == sample_name.lower():
                    return project.get('name')
    raise ValueError(f"Project was not found. project_name: [{project_name}] sample_name: [{sample_name}]")


def get_biosample_id(family_member, storage_url, auth_token, project_id):
    sample_name = family_member.get(SAMPLE_NAME, '').strip()
    sample_list = get_sample_list(storage_url, auth_token, project_id)
    for sample in sample_list:
        if sample.get("alias").lower() == sample_name.lower():
            return sample.get('name')
    return None


def get_sample_list(storage_url, auth_token, project_id):
    return get_remote_object_list(storage_url, auth_token, f'/projects/{project_id}/biosamples/')


def get_project_list(storage_url, auth_token):
    return get_remote_object_list(storage_url, auth_token, f'/projects')


def get_dataset_list_of_sample(sample_id, storage_url, auth_token):
    return get_remote_object_list(storage_url, auth_token, f'/projects/0000000000/biosamples/{sample_id}/datasets/')


def get_sequenced_files_list(sample_id, storage_url, auth_token, dataset):
    return get_remote_object_list(storage_url, auth_token,
                                  f'/projects/0000000000/biosamples/{sample_id}/datasets/{dataset.get("name")}'
                                  f'/sequenced files')


def get_sample_id_from_sample(family_member, storage_url, auth_token):
    project_list = get_project_list(storage_url, auth_token)
    # get the project id with the name as given in input
    project_id = get_project_id(family_member, storage_url, auth_token, project_list)
    # get the sample id of project id by the sample name
    return get_biosample_id(family_member, storage_url, auth_token, project_id)


def get_file_path(sample_id, dataset, file):
    return f'/projects/0000000000/biosamples/{sample_id}/datasets/{dataset.get("name")}/sequenced files/{file.get("name")}'


def get_sample_filepath_list(sample_id, dataset_list, storage_url, auth_token):
    files_paths = []
    for dataset in dataset_list:
        files = get_sequenced_files_list(sample_id, storage_url, auth_token, dataset)
        for file in files:
            file_path = get_file_path(sample_id, dataset, file)
            files_paths.append(file_path)
    return files_paths


def construct_sample_file_path_from_sample_name(family_member, host_url, storage_id, auth_token, case_group_number):
    # check if biosample cache exists
    if basespace_biosample_name_to_sequenced_files_map:
        sample_name = family_member.get(SAMPLE_NAME, '').strip()
        return basespace_biosample_name_to_sequenced_files_map.get(sample_name)

    # api storage route
    storage_url = ROUTE_GET_STORAGE_RESOURCES.format(host_url, storage_id)

    # 1) get sample id by sample name
    sample_id = get_sample_id_from_sample(family_member, storage_url, auth_token)
    if not sample_id:
        logging.error(
            f'Case group number: [{case_group_number}] could not find BioSample: {family_member.get(SAMPLE_NAME, "")}')
        return None

    # 2) get all datasets of a sample
    dataset_list = get_dataset_list_of_sample(sample_id, storage_url, auth_token)
    if not dataset_list:
        logging.error(
            f'Case group number: [{case_group_number}] could not find dataset list using BioSample id: {sample_id}')
        return None

    # 3) get list of sample file paths of all datasets
    files_paths = get_sample_filepath_list(sample_id, dataset_list, storage_url, auth_token)

    return files_paths


def get_sample_file_paths(family_member, storage_id, host_url, auth_token, case_group_number):
    samples_files_names = family_member.get(FILES_NAMES, '')
    if not is_basespace_storage_provider:
        return samples_files_names

    # basespace compatibility - check per patient
    if samples_files_names:
        # check if written with basespace id or name
        files_names = []
        for sample_file_path in samples_files_names.split(";"):
            sample_file_path = sample_file_path.strip()
            path_parts = sample_file_path.split('/')
            file_name = path_parts[-1]
            # if last file object not numeric -> convert to numeric
            if not file_name.isdigit():
                sample_file_path = convert_basespace_path_to_path_with_ids(host_url, storage_id, auth_token,
                                                                           sample_file_path.strip())
            files_names.append(sample_file_path)

        samples_files_names = files_names
    elif family_member.get(SAMPLE_NAME, ''):  # get files using sample name
        samples_files_names = construct_sample_file_path_from_sample_name(family_member, host_url, storage_id,
                                                                          auth_token, case_group_number)
    # return as in csv format: list of sample file paths, semicolon (;) separated
    return ";".join(samples_files_names)


def get_human_readable_file_name(host_url, auth_token, storage_id, sample_file_name):
    storage_url = ROUTE_GET_STORAGE_RESOURCES.format(host_url, storage_id)
    files = get_remote_object_list(storage_url, auth_token, os.path.dirname(sample_file_name.split(';')[0]))
    file_name = os.path.basename(sample_file_name.split(';')[0])
    for file in files:
        if file.get('name') == file_name:
            return file.get('alias')


def construct_sample_file_list(base_json, family_member, storage_id, host_url, auth_token, case_group_number):
    # extract sample paths
    samples_files_names = get_sample_file_paths(family_member, storage_id, host_url, auth_token,
                                                case_group_number)
    sample_files = []
    if samples_files_names:
        # set sample type for basespace only
        set_sample_type_for_basespace(auth_token, base_json, case_group_number, host_url, samples_files_names,
                                      storage_id)

        for sample_file_path in samples_files_names.split(";"):
            sample_file_path = sample_file_path.strip()
            sample_files.append({
                "filename": os.path.basename(sample_file_path),
                "path": sample_file_path,
                "size": 0,
                "storage_id": storage_id,
                "status": "uploaded",
                "vcf_column_name": family_member.get(SAMPLE_NAME, '')
            })

    bigwig = family_member.get(BIGWIG)
    if bigwig:
        sample_files.append({
            "filename": os.path.basename(bigwig.strip()),
            "path": bigwig.strip(),
            "size": 0,
            "storage_id": storage_id,
            "status": "uploaded",
            "vcf_column_name": family_member.get(SAMPLE_NAME, ''),
            "sample_type": "bw"
        })

    return sample_files


def set_sample_type_for_basespace(auth_token, base_json, case_group_number, host_url, samples_files_names, storage_id):
    if not is_basespace_storage_provider or base_json.get('sample_type'):
        return
    sample_file_name = get_human_readable_file_name(host_url, auth_token, storage_id, samples_files_names)
    set_sample_type(base_json, sample_file_name, case_group_number)


def construct_sample(family_member, sample_files, sample_type):
    sample_name = family_member.get(SAMPLE_NAME, '')
    return {
        "fastq": sample_name,
        "patient": sample_name,
        "status": "uploaded",
        "directoryPath": "",
        "sampleFiles": sample_files,
        "sample_type": "fastq" if sample_type == "fastq" else "project_vcf"
    }


def set_bam_file(family_member, sample, storage_id, host_url, auth_token):
    bam_path = family_member.get(BAM_FILE, '').strip()
    if bam_path:
        bam_location = bam_path
        if is_basespace_storage_provider:
            bam_location = get_basespace_bam_location(host_url, storage_id, auth_token, bam_path)

        sample['bam_location'] = bam_location
        sample['storage_id'] = storage_id


def get_basespace_bam_location(host_url, storage_id, auth_token, bam_path):
    path_parts = bam_path.split('/')
    file_name = path_parts[-1]
    # if last file object not numeric -> convert to numeric
    if not file_name.isdigit():
        bam_location = convert_basespace_path_to_path_with_ids(host_url, storage_id, auth_token,
                                                               bam_path.strip())
    else:
        bam_location = bam_path.strip()
    return bam_location


def validate_date_format(date):
    try:
        datetime.strptime(date, '%Y-%m-%d')
    except ValueError:
        raise ValueError(f"{date} has incorrect date format, should be YYYY-MM-DD")


def set_patient_obj(base_json, family_member, sample, host_url, auth_token):
    relationship = family_member.get(RELATION, "")
    api_relationship = "Test Subject" if relationship.lower() == 'proband' else relationship.title()

    gender = family_member.get(GENDER, '').lower()
    validate_gender(gender)
    api_gender = GENDER_DICT.get(gender[0]) if gender else 'Unknown'
    phenotypes = get_phenotype_objects(family_member, host_url, auth_token)
    api_healthy = not phenotypes

    date_of_birth = family_member.get(DATE_OF_BIRTH)

    family_member_data = {
        "fastq_sample": family_member.get(SAMPLE_NAME, ''),
        "gender": api_gender,
        "healthy": api_healthy,
        "relationship": api_relationship,
        "notes": "",
        "phenotypes": phenotypes,
        "detailed_ethnicity": get_ethnicity(),
        "id": relationship.lower()
    }

    if date_of_birth:
        validate_date_format(date_of_birth)
        family_member_data['date_of_birth'] = date_of_birth

    if relationship.lower() == 'sibling':
        base_json['patients']['other'].append(family_member_data)
    else:
        base_json['patients'][relationship] = family_member_data

    base_json['samples'].append(sample)


def validate_gender(gender):
    if not gender.lower() in GENDER_DICT.keys():
        raise ValueError(f"Gender must be specified: M|F|U - male, female or unknown. Current value: [{gender}]")


def validate_relationship(relationship, existing_family_members):
    if relationship.lower() not in RELATIONS:
        raise ValueError(f"Relation is not valid. value: [{relationship}]")
    if relationship != 'sibling' and relationship in existing_family_members:
        raise ValueError(f"{relationship} must be unique. There is more than one [{relationship}] in this family.")
    else:
        existing_family_members.add(relationship)


def set_boost_genes(family_member, base_json):
    boost_genes = family_member.get(BOOST_GENES, '')
    base_json['boostGenes'] = True if boost_genes.lower() == 'true' else False


def set_gene_list(family_member, base_json):
    gene_list_id = family_member.get(GENE_LIST_ID)
    if gene_list_id:
        try:
            gene_list_id = int(gene_list_id)
        except ValueError:
            raise ValueError(f"Gene list id is not a valid number. value: {gene_list_id}")

        base_json['gene_list']['type'] = 'existing'
        base_json['gene_list']['id'] = gene_list_id
        base_json['gene_list']['visible'] = True


def set_kit(family_member, base_json):
    kit_id = family_member.get(KIT_ID)
    if kit_id:
        try:
            kit_id_int = int(kit_id)
            base_json['sequence_info'] = {
                'kit_id': kit_id_int
            }
        except ValueError:
            raise ValueError(f"Kit id is not a valid number. value: {kit_id}")


def set_case_type(base_json, family_member):
    case_type = family_member.get(CASE_TYPE)
    if case_type:
        base_json['type'] = case_type


def does_proband_have_siblings(base_json):
    if base_json['patients']['other']:
        for s in base_json['patients']['other']:
            if s.get('relationship').lower() == 'sibling':
                return True
    return False


def get_all_siblings_idxs(base_json):
    siblings = []
    for i in range(len(base_json['patients'].get('other'))):
        if base_json['patients']['other'][i].get('relationship').lower() == 'sibling':
            siblings.append(i)
    return siblings


def set_parents_id(base_json):
    mother_id = base_json['patients'].get('mother', {}).get('id', '')
    proband_has_sibling = does_proband_have_siblings(base_json)
    siblings = get_all_siblings_idxs(base_json)
    if mother_id:
        base_json['patients']['proband']['mother'] = mother_id
        if proband_has_sibling:
            for s in siblings:
                base_json['patients']['other'][s]['mother'] = mother_id
    father_id = base_json['patients'].get('father', {}).get('id', '')
    if father_id:
        base_json['patients']['proband']['father'] = father_id
        if proband_has_sibling:
            for s in siblings:
                base_json['patients']['other'][s]['father'] = father_id


def set_selected_preset(family_member, base_json):
    selected_preset = family_member.get(SELECTED_PRESET)
    if selected_preset:
        base_json['selected_preset_set'] = selected_preset


def set_due_date(family_member, base_json):
    due_date = family_member.get(DUE_DATE)
    if due_date:
        validate_date_format(due_date)
        base_json['due_date'] = due_date


def set_notes(family_member, base_json):
    notes = family_member.get(CLINICAL_NOTES)
    if notes:
        base_json['notes'] = notes


def set_label(family_member, base_json):
    labels = family_member.get(LABELS)
    if labels:
        try:
            labels_list = list(set(int(i.strip()) for i in labels.split(';')))
            base_json['labels'] = labels_list
        except ValueError:
            logging.exception("One of the labels is not a valid number, creating a case without labels")


def get_case_data(family_members, storage_id, host_url, auth_token, case_group_number):
    logging.info(f"Case group number: [{case_group_number}] Start generating ANC payload json")

    # get empty/default ANC payload json
    base_json = get_base_json()
    should_upload = False
    existing_family_members = set()

    # traverse family members
    for family_member in family_members:
        relationship = family_member.get(RELATION, "").lower()
        validate_relationship(relationship, existing_family_members)
        # case settings - only by proband record
        if relationship == 'proband':
            should_upload = set_proband_object(base_json, case_group_number, family_member)

        # set patient's sample/files
        sample_files = construct_sample_file_list(base_json, family_member, storage_id, host_url, auth_token,
                                                  case_group_number)
        sample = construct_sample(family_member, sample_files, base_json.get('sample_type'))
        set_bam_file(family_member, sample, storage_id, host_url, auth_token)
        set_patient_obj(base_json, family_member, sample, host_url, auth_token)

    set_parents_id(base_json)
    enrich_disease_phenotypes = False
    include_incidental_parents = False

    logging.info(f"Case group number: [{case_group_number}] payload json is created")

    return base_json, should_upload, enrich_disease_phenotypes, include_incidental_parents


def set_proband_object(base_json, case_group_number, family_member):
    # set should upload indicator
    execute_now = family_member.get(EXECUTE_NOW, 'true').lower()
    if execute_now not in ['true', 'false', '']:
        raise ValueError(f"Execute now is not valid. value: [{execute_now}]")
    should_upload = not (execute_now in ['', 'true'])

    # set samples and case types
    set_sample_type(base_json, family_member.get(FILES_NAMES, ''), case_group_number)
    set_case_type(base_json, family_member)
    set_boost_genes(family_member, base_json)
    set_gene_list(family_member, base_json)
    set_kit(family_member, base_json)
    set_selected_preset(family_member, base_json)
    set_due_date(family_member, base_json)
    set_notes(family_member, base_json)
    set_label(family_member, base_json)
    return should_upload


def create_case(test_data, auth_token, host_url, case_group_number):
    logging.info(f'Case group number: [{case_group_number}] - before sending ANC payload to api server')
    url = ROUTE_CREATE_CASE.format(host_url)
    res = requests.post(url=url, json=test_data, headers={'Authorization': auth_token}, timeout=600)

    # case was not created
    if res.status_code != 201:
        logging.error(
            f'Case group number: [{case_group_number}] - Case creation failed. status code: {str(res.status_code)} '
            f'server response: {res.text}  payload sent: {test_data}')
        return False

    # case was created
    response = res.json()
    logging.info(f'Case group number: [{case_group_number}] - Case creation succeeded. server response: {response}')
    return True


def is_json(info_file_path):
    return pathlib.Path(info_file_path).suffix.lower() == ".json"


def construct_case_data(info_file_path, storage_id, host_url, auth_token):
    counter = 0
    with open(info_file_path, mode='r') as info_file:
        info_reader = json.load(info_file) if is_json(info_file_path) else csv.DictReader(info_file)
        for case_group_number, info_lines in groupby(info_reader, key=lambda x: x[CASE_GROUP_NUMBER]):
            counter += 1
            if not case_group_number:
                logging.error(f'Case group number was not set. continue.')
                yield None, case_group_number
                continue

            logging.info(f'Case group number: [{case_group_number}] - generating test data')
            try:
                case_data, should_upload, enrich_disease_phenotypes, include_incidental_parents = \
                    get_case_data(info_lines, storage_id, host_url, auth_token, case_group_number)
            except Exception as e:
                logging.exception(f"Case group number: [{case_group_number}] - Error: {e}. Continue to next case.")
                yield None, case_group_number
                continue

            # ANC payload
            test_data = {
                "enrich_disease_phenotypes": enrich_disease_phenotypes,
                "include_incidental_parents": include_incidental_parents,
                "should_upload": should_upload,
                "test_data": case_data
            }
            yield test_data, case_group_number

        logging.info(f'cases count: {counter}')


def validate_organization(api_auth_token, host_url, user_name):
    url = ROUTE_USERS.format(host_url)
    res = requests.get(url=url, headers={'Authorization': api_auth_token})
    if res.status_code != 200:
        logging.error(f'Status_code: {res.status_code} text: {res.text}')
        res.raise_for_status()
    else:
        response = res.json()
        for user in response['hits']:
            if user['email'] != user_name:
                continue
            org_id = user.get('organization_id')
            loop = True
            while loop:
                logging.info(f"Current user is on organization id: {org_id}")
                answer = input("Is it OK to continue? Enter 'yes' or 'no'")
                if answer.lower() == "yes":
                    loop = False
                elif answer.lower() == "no":
                    exit(f"User: {user_name} org id: {org_id} - user wishes to exit script")
                else:
                    print("Please enter yes or no.")
            if not loop:
                break


def create_cases(api_auth_token, info_file_path, storage_id, host_url):
    new_cases_counter = 0
    failed_cases_counter = 0
    case_index = 0

    for family_test_data, case_group_number in construct_case_data(info_file_path, storage_id, host_url,
                                                                   api_auth_token):
        case_index += 1
        if not family_test_data:
            failed_cases_counter += 1
            logging.exception(f"Case group number: [{case_group_number}]. Error in creating case data. "
                              f"Continue to next case")
            continue
        try:
            case_created = create_case(family_test_data, api_auth_token, host_url, case_group_number)
            if case_created:
                new_cases_counter += 1
            else:
                failed_cases_counter += 1
        except Exception as e:
            logging.exception(f"Case group number: [{case_group_number}]. Error in creating new case: {e}. "
                              f"Continue to next case.")
            failed_cases_counter += 1

    logging.info(f"create_cases is done. total cases: {case_index} "
                 f"new_cases_counter: {new_cases_counter} "
                 f"failed_cases_counter {failed_cases_counter}")


def generate_extractor(args):
    return PhenotypeExecutor(args.phenotype_extractor_type,
                             args.phenotypes_domain,
                             args.phenotypes_username,
                             args.phenotypes_password,
                             args.phenotypes_secret)


def prepare_biosample_cache(host_url, storage_id, auth_token):
    global basespace_biosample_name_to_sequenced_files_map

    # api storage route
    storage_url = ROUTE_GET_STORAGE_RESOURCES.format(host_url, storage_id)

    # 1) get all projects
    project_list = get_project_list(storage_url, auth_token)

    # 2) loop over all projects
    for project in project_list:
        # for each project get the list of samples
        sample_list = get_sample_list(storage_url, auth_token, project.get('name'))
        for sample in sample_list:
            # for each sample get the dataset list and then get the list of sample files
            sample_id = sample.get('name')
            dataset_list = get_dataset_list_of_sample(sample_id, storage_url, auth_token)
            samples_files_names = get_sample_filepath_list(sample_id, dataset_list, storage_url, auth_token)
            sample_name = sample.get('alias')
            # map the list of sample files to the sample name
            basespace_biosample_name_to_sequenced_files_map[sample_name] = samples_files_names


def main():
    global is_basespace_storage_provider, phenotype_extractor_type, phenotype_extractor

    # read command line arguments
    args = parse_args()

    if args.version:
        print(f"version: {VERSION}")
        return

    # logger
    setup_logger(log_folder='create_batch_cases_v2_logs', log_prefix='create_batch_cases_v2')

    # storage and phenotype provider
    is_basespace_storage_provider = args.is_basespace
    phenotype_extractor_type = args.phenotype_extractor_type
    if phenotype_extractor_type:
        phenotype_extractor = generate_extractor(args)

    logging.info(f"Create_batch_script_v2 version: {VERSION} - "
                 f"Reads patients' metadata form an input csv/json file and sends a "
                 f"create case request, per patient-group.")

    # get the api token
    api_auth_token = get_authentication_token(args.username, args.password, args.host_url)

    if DEBUG_VERSION:
        validate_organization(api_auth_token, args.host_url, args.username)

    # prepare the basespace cache
    if is_basespace_storage_provider:
        prepare_biosample_cache(args.host_url, args.storage_id, api_auth_token)

    # parse the input file and generate the api requests for case creation
    create_cases(api_auth_token, args.info_file_path, args.storage_id, args.host_url)

    logging.info("create_batch_cases_v2 is finished")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

"""
Utilities for the GAPP (was Rapidomics) project, including tools to interact 
with Phenotips, Illumina's TruSight Software Suite (TSS), BaseSpace Sequence
Hub (BSSH).
"""

import sys, os, warnings
import logging
import re
import json
import requests
import pandas as pd

__version__ = 0.1

class Configurator:
    """
    Parse configuration files.
    """
    def __init__(self, config_file=os.path.expanduser("~/.illumina/gapp_conf.json")):
        """
        Load settings from config_file (JSON). Uses default if none provided. 
        Returns self, where attributes are settings parsed from config_file.
        """
        self.file = config_file
        try:
            with open(self.file, 'r') as cfg:
                self.configs =json.load(cfg)
        except FileNotFoundError as errmsg:
            sys.exit(errmsg)
        except Exception as errmsg:
            sys.exit(errmsg)
        else:
            self.instance         = self.configs['instance']
            self.domain           = self.configs['X-ILMN-Domain']
            self.workgroup        = self.configs['X-ILMN-Workgroup']
            self.token            = self.configs['X-Auth-Token']
            self.testDefinitionId = self.configs['testDefinitionId']
            
            # Authorizations for BaseSpace Sequence Hub (BSSH) are written by
            # Illumina's CLI tool `bs` to file ~/.basespace/*.cfg
            #
            self.bs_server = self.configs['bs_apiServer']
            self.bs_token  = self.configs['bs_accessToken']

            # Configurations for Phenotips
            #
            self.phenotips_server = self.configs['X-Gene42-Server']
            self.phenotips_auth   = self.configs['X-Gene42-Auth']
            self.phenotips_secret = self.configs['X-Gene42-Secret']

            # Configurations for REDCap
            #
            self.redcap_server = self.configs['REDCap-Server']
            self.redcap_token  = self.configs['REDCap-Token']

            # Configurations for Emedgene
            #
            self.emg_username    = self.configs['EMG-Username']
            self.emg_password    = self.configs['EMG-Password']
            self.emg_prag_server = self.configs['EMG-PRAG-Server']
            self.emg_eval_server = self.configs['EMG-EVAL-Server']


class REDCap:
    """
    Return object for [REDCap](https://tacc-redcap.bic.mni.mcgill.ca/) API
    """
    def __init__(self, config_file=os.path.expanduser("~/.illumina/gapp_conf.json")):
        """
        Load settings from config_file, if provided. Define instance vars to
        provide more readable access to settings in dict "configs".
        """
        configs     = Configurator()
        self.server = configs.redcap_server
        self.token  = configs.redcap_token

    def get_record_id(self, q1k_id):
        """
        Get REDCap record_id for patient by `q1k_id`.
        - q1k_id: [str] of the following format e.g.: Q1K_HSJ_10050_P.
        - Returns [str] e.g.: '50', or None (error).
        """
        filterLogic = '[q1k_proband_id_1]="' + q1k_id + '" or [q1k_relative_idgenerated_1]="' + q1k_id + '"'
        data = {
            'token': self.token,
            'content': 'record',
            'action': 'export',
            'format': 'json',
            'type': 'flat',
            'csvDelimiter': '',
            'fields[0]': 'record_id',
            'rawOrLabel': 'raw',
            'rawOrLabelHeaders': 'raw',
            'exportCheckboxLabel': 'false',
            'exportSurveyFields': 'false',
            'exportDataAccessGroups': 'false',
            'returnFormat': 'json',
            'filterLogic': filterLogic
        }
        response = requests.post(self.server, data=data)
        if response.status_code == 200:
            return response.json()[0]['record_id']
        else:
            logging.error('HTTP Status: ' + str(response.status_code))
            logging.debug(json.dumps(response.json(), indent=2))
            return None

    def get_hpo(self, q1k_id):
        """
        Get HPO terms for `q1k_id`.
        - q1k_id: [str] of the following format e.g.: Q1K_HSJ_10050_P.
        - Returns [list] e.g.: '50', or None (error).
        """
        hpos = []
        record_id = self.get_record_id(q1k_id)
        data = {
            'token': self.token,
            'content': 'record', # required
            'action': 'export',  # default? Not required
            'format': 'json',    # required
            'rawOrLabel': 'raw', # 'label' will return HPO Label instead of HPO ID
            'events[0]': 'direct_assessment_arm_1', # Get all events if not provided
            'records[0]': record_id #,
            # 'fields[n]' NOT specified seem to return all fields
            # 'fields[0]': 'record_id', # to include record_id in the json output
            # 'fields[1]': 'ghf_vis_hpo_1', (...) 'fields[9]': 'ghf_dely_gmd_hpo'
        }
        response = requests.post('https://tacc-redcap.bic.mni.mcgill.ca/api/',data=data)
        if response.status_code == 200:
            for value in response.json()[0].values():
                hpos.append(value) if value.startswith('HP:') else None
        else:
            logging.error('REDCap.get_hpo() HTTP Status: ' + str(response.status_code))
            logging.debug(json.dumps(response.json(), indent=2))
        return ';'.join(hpos)


class Phenotips:
    """
    Return an object, to interact with Phenotips.
    """
    def __init__(self, config_file=os.path.expanduser("~/.illumina/gapp_conf.json")):
        """
        Load settings from config_file, if provided. Define instance vars to
        provide more readable access to settings in dict "configs".
        """
        configs      = Configurator()
        self.server  = configs.phenotips_server
        self.auth    = configs.phenotips_auth
        self.secret  = configs.phenotips_secret
        self.headers = {
            'Authorization'  : configs.phenotips_auth,
            'X-Gene42-Secret': configs.phenotips_secret
        }
    

    def get_patient(self, pid):
        """
        Get patient information by Phenotips ID (`pid`).
        - `pid`: `str` of the following format e.g.: P0000001).
        - Returns A dict (requests.json() object) for `pid`, or None (error).
        """
        try:
            if re.match(r"P\d{7}", pid):
                url = self.server + '/rest/patients/' + pid
            else:
                print(f"Wrong format for PID: '{pid}'. Should be like 'P0000001'.")
                return(None)
        except TypeError:
            print(f"Phenotips ID (PID) '{pid}' is not of type 'str'")
            return(None)
        except Exception as errmsg:
            print(errmsg)
            return(None)

        response = requests.get(url, headers=self.headers)
        if response.status_code != 200:
            return(None)
        else:
            return(response.json())


    def get_patient_by_eid(self, eid):
        """
        Get patient by External ID (`eid`).
        - `eid`: `str` for the External IDentifier, e.g.: 72111017.
                 `eid` is actually a truncated RAMQ number (without letters)
        - Returns A dict (requests.json() object) for `eid`, or None (error).
        """
        url = self.server + '/rest/patients/eid/' + eid
        response = requests.get(url, headers=self.headers)
        if response.status_code != 200:
            response.raise_for_status()
            return(None)
        else:
            return(response.json())


    def get_patient_by_mrn(self, mrn):
        """
        Get patient by a labeled external ID (Medical Record Number, MRN).
        https://docs.phenotips.com/reference/getpatientbylabeledeid
        - `mrn`: `str` for the labeled external IDentifier
        - Returns Requests.json() object for `mrn` or None (error)
        """
        label = 'MRN'
        url   = self.server + f"/rest/patients/labeled-eid/{label}/{mrn}"
        try:
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
        except requests.exceptions.ConnectionError as err:
            raise SystemExit(err)
        except requests.exceptions.HTTPError as err:
            return(None)
        else:
            return(response.json())


    def parse_hpo(self, phenotips_json):
        """
        Retrieve HPO terms observed in `phenotips_json`.
        - `phenotips_json`: Phenotips patient data (dict for JSON).
        - Returns: List of dicts for HPO terms that are 'observed'='yes',
          or an empty list if something went wrong. E.g.:
          [{'id': 'HP:0001561', 'label': 'Polyhydramnios'}, {}, {},...] or [].
        """
        hpos = []

        # List of HPO terms in list of dicts `phenotips_json['features']`, eg:
        # [{'id': 'HP:0001250', 'label': 'Seizures', 'type': 'phenotype', 'observed': 'yes'},
        #  {'id': 'HP:0001251', 'label': 'Ataxia',   'type': 'phenotype', 'observed': 'no'}, 
        #  {etc...}]
        # Report only HPO features that are 'observed', in a list that can be
        # used to create Cases for TSS (eg: [{'code': 'HP:0000589', 'source': 'HPO'}, {...}]
        #
        for feature in phenotips_json['features']:
            if feature['observed'] == 'yes':
                hpos.append({'id': feature['id'], 'label': feature['label']})
        return(hpos)


    def get_hpo(self, pid):
        """
        Get HPO terms observed for this Phenotips ID (`pid`, must be a `str`).
        - `pid`: `str` of the following format e.g.: P0000001).
        - Returns: List of dicts for HPO terms that are 'observed'='yes'
          Returns an empty list if something went wrong. E.g.:
          [{'id': 'HP:0001561', 'label': 'Polyhydramnios'}, {}, {},...] or [].
        """
        patient = self.get_patient(pid)
        return(self.parse_hpo(patient))


    def get_hpo_old_from_db_folder(self, pid, db='/staging2/data/Illumina/TSS/2021-03-08/json/pheno_json'):
        """
        Get HPO terms observed for this Phenotips ID (`pid`, must be a `str`),
        from the backup of the old Phenotips system.
        - `pid`: `str` of the following format e.g.: P0000001).
        - `db` : Old database backup, a folder containing all the patients
          record files (as JSON). Default location is 
          spxp-app02://staging2/data/Illumina/TSS/2021-03-08/json/pheno_json/
        - Returns: List of dicts for HPO terms that are 'observed'='yes'
          Returns an empty list if none found, or if something went wrong.
        """
        # Collect HPO terms that are marked 'yes' for 'observed'.
        #
        hpos = []
        
        # Use os.path.isfile() to find the json file to `pid` in `db`.
        # File-naming convention under the `db` folder is "phenotipsid.json",
        # e.g. './pheno_json/P0000808.json'.
        #
        filename = db + os.sep + pid + '.json'
        if os.path.isfile(filename):
            with open(filename, 'rb') as fh:
                patient = json.load(fh)
            # Parse the json to build list of HPO terms
            #
            for feature in patient['features']:
                if feature['observed'] == 'yes':
                    hpos.append({'id': feature['id'], 'label': feature['label']})
        else:
            warnings.warn(f"WARNING: '{filename}' is not a file {filename}")
        return(hpos)


    def load_old_pids_to_dict(self, db='/staging2/data/Illumina/TSS/2021-03-08/json/pheno_json', outfile=None):
        """
        Load as dict all the old PID entries from `db` backup. Do this once.
        - `db` : [str] Old database backup folder containing all the patients
          record files (as JSON). Default location is 
          spxp-app02://staging2/data/Illumina/TSS/2021-03-08/json/pheno_json/
        - `outfile`: [str] Filename to JSON output. Default=None (not saved).
          Backup as 'pheno_json-extracted_hpos.json'
        - Returns: a dict for HPO terms that are 'observed'='yes'
          {"PID": {'id': [hpo1,hpo2,...], 'label': [label1,label2,...], },...}
        """
        pids = {}
        pid_files = os.listdir(db)
        for pid_file in pid_files:
            pid = os.path.splitext(pid_file)[0]
            pids[pid] = self.get_hpo_old_from_db_folder(pid, db=db)

        if outfile is not None:
            with open(outfile, 'w') as fh:
                fh.write(json.dumps(pids, indent=2))

        return(pids)


    def get_hpo_old(self, pid, db='/staging2/data/Illumina/TSS/2021-03-08/json/pheno_json.json'):
        """
        Get HPO terms for PID from the file 'old_phenotips_extracted_hpos.json'
        which contains the old data extracted from the backup using 
        `self.load_old_pids_to_dict(
            db='/staging2/data/Illumina/TSS/2021-03-08/json/pheno_json', 
            outfile='old_phenotips_extracted_hpos.json')`
        - `pid`: `str` of the following format e.g.: P0000001).
        - `db` : JSON file with PID and HPOs extracted from the old Phenotips. 
          spxp-app02://staging2/data/Illumina/TSS/2021-03-08/json/pheno_json/
        - Returns: List of HPO terms that are 'observed'='yes' for PID.
          Returns an empty list if none found, or if something went wrong.
        """
        with open(db) as fh:
            pids = json.loads(fh.read())
        return(pids[pid])


class TSS:
    """
    Return a TSS object, with API functionalities.
    """
    def __init__(self, config_file=os.path.expanduser("~/.illumina/gapp_conf.json")):
        """
        Load settings from config_file, if provided. Define instance vars to
        provide more readable access to settings in dict "configs".
        """
        configs               = Configurator(config_file)
        self.instance         = configs.instance
        self.domain           = configs.domain
        self.workgroup        = configs.workgroup
        self.token            = configs.token
        self.testDefinitionId = configs.testDefinitionId
        self.server           = 'https://' + self.instance
        self.headers = {
            'X-Auth-Token'    : configs.token,
            'X-ILMN-Domain'   : configs.domain,
            'X-ILMN-Workgroup': configs.workgroup
        }
        self.cases = pd.DataFrame()
    
    def list_cases(self):
        """
        GET a list of cases, with correspondances between 'id' and 'CaseID'
        ('displayId').
        Sets list of correspondances as a pandas DataFrame object in attribute
        `self.cases` (a dict) and returns it.
        """
        endpoint = '/crs/api/v2/cases/search?displayId=%25'
        url = self.server + endpoint
        response = requests.get(url, headers=self.headers)
        
        # response is a JSON string. Use json.loads() to convert into dict of
        # 11 keys. Key 'content' holds a list of 'id' to 'displayId' (CaseID).
        # [{'id': '01d4d79b-a33c...', 'displayId': 'GM...'}, {}, {},...]
        #
        json_data = json.loads(response.text)
        self.cases = pd.json_normalize(json_data['content'])
        return(json_data['totalElements'])
    
    def search_case(self, displayId):
        """
        GET id for a CaseID
        - displayId: Search string. Allows partial search 
          e.g. ILM-ABC-234, ABC, 234 will include ILM-ABC-234 in result
        - Returns: A list of `CaseID` matching `displayId`.
        """
        endpoint = '/crs/api/v2/cases/search?displayId=' + displayId
        url = self.server + endpoint
        response = requests.get(url, headers=self.headers)
        jsondata = json.loads(response.text)
        case_ids = []
        for case in jsondata['content']:
            case_ids.append(case['id'])
        return(case_ids)

    def get_case(self, displayId):
        """
        GET the JSON for a CaseID, using the displayID
        - displayId: String, _e.g._ GM200578
        - Returns: JSON for `displayId` (as a dict). Ex: 
            - {'id': 'None found!', 'displayId': 'HM210196'}
            - {'id': "More than one found ['60eb7387-585b-4b76-b072-433f7359cdd2', 'b25da3c5-ed37-43de-aff3-1b796f2b1ce5', 'b3eeeb07-945e-469d-8846-b7d01d29011a']", 'displayId': 'HM210011'}
            - {'id': 'cf2136dd-1b32-4a55-99d2-af941cf355dc', 'displayId': 'HM210221', 'createdDate': 'blah', ...}
        """
        caseid = self.search_case(displayId)
        if len(caseid) == 0:
            #return(f"{{'id':'None found!', 'displayId':'{displayId}'}}")
            return {'id': 'None found!', 'displayId': displayId}
        elif len(caseid) > 1:
            return  {'id': f'More than one found {caseid}', 'displayId': displayId}
        elif len(caseid) == 1:
            endpoint = '/crs/api/v1/cases/' + caseid[0]
            url = self.server + endpoint
            response = requests.get(url, headers=self.headers)
            jsondata = json.loads(response.text)
            return(jsondata)


    def delete_case(self, displayId):
        """
        Search and DELETE Case by displayId. Deletes case if only one found.
        - displayId: Case to delete 
        - Returns: ['status', [list, of, ids, found]] where 'status' can be:
            - Deleted: OK
            - Failed: Raise warning, one found but could not delete from TSS.
            - Absent: displayId not found in TSS
            - Multiple: More than one corresponding displayId found
        """
        results = self.search_case(displayId)
        if len(results) == 1:
            endpoint = '/crs/api/v1/cases/' + results[0]
            url = self.server + endpoint
            response = requests.request("DELETE", url, headers=self.headers)
            if response.ok:
                status = 'Deleted'
            else:
                status = 'Failed'
                warnings.warn(f'Failed to delete on TSS: {response.status_code}')
        elif len(results) == 0:
            status = 'Absent'
        else:
            status = 'Multiple'
        return([status, results])

    def create_case(self, case):
        """
        POST `case`, a JSON string that describes a Case, to the TSS web API.
        Returns the API's response, as a 'Requests' object.
        """
        endpoint = '/crs/api/v1/cases/?forceOverwrite=true'
        url      = self.server + endpoint
        headers  = self.headers
        headers['Content-Type'] = 'application/json'
        
        # Case has to be submitted as a JSON str, else TSS will return an HTTP
        # error code 400. Cases as dict are serialized into JSON str.
        #
        if type(case) is str:
            response = requests.post(url, headers=headers, data=case)
        elif type(case) is dict:
            try: 
                response = requests.post(url, headers=headers, data=json.dumps(case))
            except TypeError as err:
                print("Convert dict to json failed. Are you sure you passed a dict?")
                print(err)
            except Exception as err:
                print(f"Something unexpected happened. {err}")
        else:
            raise TypeError(f"Expected case as a JSON string, but got a {type(case)}")
        
        return(response)
    
    def download_files(self, displayId, pattern='vcf.gz', matchExactPath='true'):
        """
        Download 'file' type from TSS for a given displayId (case name).
        - displayId: Case for which we wish to download the file
        - pattern: Pattern in filename. e.g. 'vcf.gz', 'vcf.gz.tbi';
        Returns a list of tuples [(gds, url), (gds, url), ...]. 
          Use url to download file at gds in TSS
        """
        endpoint = '/crs/api/v1/files?'
        headers  = self.headers
        headers['Content-Type'] = 'application/json'
        paths = [] # to hold (gds, url) tuples 

        # A full GDS path is required to download files for a case.
        # Find the TSS internal 'id' for a case ('displayId') to get the case
        # definition, as JSON. The full gds:// is deeply nested in the JSON.
        #
        case = self.get_case(displayId)

        # Within 'case', 'ingestionResult' contains lists of all gds paths as
        # JSON but is encoded as a string. Decode to parse.
        #
        try:
            ingestionResult = json.loads(case['ingestionResult'])
        except TypeError as err:
            print(f"ERROR while retrieving case[ingestionResult] for {displayId}:\n\t{err}")
        else:
            # Within ingestionResult, outputFileTables is a list of result types,
            # stored as dicts, where 'name' identifies the positional index: 
            # [0:Alignments, 1:Variants, 2:Metrics, 3:Visualizations, 4:Reports]
            # Each result type contain a list of files, and each file is a dict
            # describing the metadata: {type, path, sizeInBytes, tags}.
            #
            paths = [] # to hold (gds, url) tuples 
            for table in ingestionResult['result']['outputFileTables']:
                for file in table['files']:
                    if file['path'].find(pattern) != -1:
                        url = self.server + endpoint + f"includePresignedUrl=true" + f"&matchExactPath={matchExactPath}" + f"&path={file['path']}"
                        paths.append((file['path'], url))
            #response = requests.request("GET", url, headers=headers)
        return(paths)

    def connect_biosamples_to_case(self):
        """
        TODO: Do we really need to use the `bs` CLI, or is there a web API?
        """
        pass


class BSSH:
    """
    Return a BSSH object for interacting with Illumina BaseSpace Sequence Hub,
    """
    def __init__(self, config_file=os.path.expanduser("~/.illumina/gapp_conf.json")):
        """
        Load settings from config_file, if provided. Define instance vars to
        provide more readable access to settings in dict "configs".
        """
        configs      = Configurator(config_file)
        self.server  = configs.bs_server
        self.token   = configs.bs_token
        self.headers = {'Authorization': f'Bearer {configs.bs_token}'}
        # self.headers = {'x-access-token': f'{token}'} # Also works


    def get_biosampleid(self, biosamplename):
        """
        Get BSSH ID for biosamplename
        - `biosamplename`: `str`, name of biosample (e.g. '28849')
        - Returns id as `int` or NoneType if not found 
        """
        endpoint = '/v2/biosamples/'
        url      = self.server + endpoint
        payload  = {'biosamplename': f"{biosamplename}"}
        response = requests.get(url, headers=self.headers, params=payload)
        response.raise_for_status
        # TODO: Warn if response.json().get('Paging')["TotalCount"] != 1
        return response.json().get('Items')[0]['Id']


    def get_datasets(self, biosampleid):
        """
        Get BSSH dataset ID, project ID and project Name for `biosampleid`.
        - `biosampleid`: [str], Id of biosample
        - Returns: [list] of tuples ['ds.id', 'project.id', project.name]
        """
        endpoint = '/v2/datasets/'
        url      = self.server + endpoint
        payload  = {'inputbiosamples': {biosampleid}} # , 'datasettypes': ''}
        # FastQ uploaded using CLI has DatasetTypes.ID 'illumina.fastq.v1.8' 
        # while the ones created by BCL Convert have the type 'common.fastq'.

        response = requests.get(url, headers=self.headers, params=payload)
        response.raise_for_status

        items  = response.json().get('Items')
        counts = response.json().get('Paging')['TotalCount']

        datasets = []
        if len(items) != counts:
            warnings.warn(f"WARNING: Found {len(items)} datasets but expected {counts} for {biosampleid}")
        for item in items:
            datasets.append((item['Id'], item['Project']['Id'], item['Project']['Name']))
        if len(datasets) != 1:
            warnings.warn(f"WARNING: Found more than one dataset for {biosampleid}")
        return datasets


    def get_sequenced_files(self, biosample):
        """
        Get paths to sequenced files (FASTQ paths) for biosample
        - `biosample`: biosamplename [str]
        - Returns: List of BSSH paths to sequenced files (FASTQ) [list of str]
        For each sequenced file, we get a path that looks like this:
        /projects/###/biosamples/###/datasets/###/sequenced files/###
        """
        fastqs = []
        biosampleid = self.get_biosampleid(biosample)
        datasets    = self.get_datasets(biosampleid)
        # projectid   = '0000000000'
        # '0000000000' works but Illumina recommends using ProjectId identifier
        # which can be found in the JSON returned for ds.dataseid 
        #
        # For each biosampleid, there's one or more datasets (runs or lanes?)
        # wich, in turn, has one or more files (sequencing lanes and reads)
        # https://api.cac1.sh.basespace.illumina.com/v2/datasets/ds.16ef9da98da745bf9f51c57e89c282f6/files/
        #
        for dataset in datasets:
            datasetid, projectid, projectname = dataset
            endpoint = f"/v2/datasets/{datasetid}/files"
            url      = self.server + endpoint
            payload  = {'limit': 100}
            response = requests.get(url, headers=self.headers, params=payload)
            response.raise_for_status
            for item in response.json().get('Items'):
                #print(json.dumps(item, indent=2))
                # Name of "biosample" works but EMG recommends "biosampleid"
                # fastq = f"/projects/0000000000/biosamples/{biosample}/datasets/{datasetid}/sequenced files/{item['Id']}"
                fastq = f"/projects/{projectid}/biosamples/{biosampleid}/datasets/{datasetid}/sequenced files/{item['Id']}"
                fastqs.append(fastq)

        return fastqs 


class Emedgene:
    """
    Return object for interacting with Emedgene's API
    """
    def __init__(self, config_file=os.path.expanduser("~/.illumina/gapp_conf.json")):
        """
        Load settings from config_file, if provided. Define instance vars to
        provide more readable access to settings in dict "configs".
        """
        configs          = Configurator(config_file)
        self.username    = configs.emg_username
        self.password    = configs.emg_password
        self.prag_server = configs.emg_prag_server
        self.eval_server = configs.emg_eval_server


    def authenticate(self):
        """
        Returns an authorization token.
        N.B. The Authorization header expires after 8H, after that, requests 
        will return an error code 403. To resolve, re-do the Login procedure to
        get a new token.
        """
        # TODO: Add different domain servers
        url      = f"{self.prag_server}/api/auth/api_login/"
        payload  = f'{{"username": "{self.username}", "password": "{self.password}"}}'
        headers  = {'Content-Type': 'application/json'}
        response = requests.request("POST", url, headers=headers, data=payload)
        return response.json()["Authorization"]


    def get_emg_id(self, sample):
        """
        Returns EMG identifier for Sample
        - `sample`: Sample Names (ex.: GM232823, 24-00666-T1, MO-24-003708,...)
        - Returns : [str] ex.: EMG107903188, None (not found) or HTTPErrorCode
        """
        # TODO: Add different domain servers
        url = f"{self.prag_server}/api/sample/?query={sample}&sampleType=fastq"
        resp = requests.get(url, headers={'Authorization': self.authenticate()})
        if resp.status_code == 200:
            if resp.json()['total'] == 1:
                return resp.json()['hits'][0]['note']
            elif resp.json()['total'] == 0:
                return None
            else:
                logging.warning(f"More than one result found: {resp.json()['total']}")
                return resp.json()
        elif resp.status_code == 401 or resp.status_code == 403:
            logging.warning(f"Unauthorized: please authenticate yourself")
            return resp.status_code
        else:
            logging.warning(f"While fetching EMG ID, got the HTTP Error Code: [{resp.status_code}]\n{resp.text}")
            return resp.status_code


    def submit_emg_case(self, case_json, server):
        """
        Submit to Emedgene a JSON that describes a Case for creation. 
        - case_json: [str] path/to/case_json
        - Return:    [obj] Response object from a `requests()` POST.
        """
        url = f"{server}/api/cases/v2/cases/"
        with open(case_json, 'r') as fh:
            payload = json.load(fh)
            resp = requests.post(url, headers={'Authorization': self.authenticate()}, json=payload)
        return resp.json()


def test():
    """
    Quick and dirty testing
    """
    pho = Phenotips()
    pho.get_hpo_old('P0000808')


def main():
    print("\nDone.\n")


if __name__ == '__main__':
    test()
    main()

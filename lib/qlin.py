#!/usr/bin/env python3

"""
Utilities for the QLIN API that creates/stores analyses, patients and sequencings in the system

Import as a package from other modules:

>>> python
from qlin import qlin
q = qlin(url)

As a standalone app, get the standard usage:

>>> shell
clin.py --help

"""

import requests
import json
import sys
import os
import argparse
import re
import pandas as pd
import random
import string
import logging
#To format dates properly
from datetime import datetime
#importing necessary functions from dotenv library
from dotenv import load_dotenv, dotenv_values

# Imports from current lib
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from lib.nanuq import Nanuq

__version__ = "0.3"

# Get a logger for this module using __name__ for hierarchical logging
logger = logging.getLogger(__name__)

class ArgumentException (Exception):
    """
    Defines an ArgumentException
    """
    def __init__(self, message):
        super().__init__(message)

class APIException (Exception):
    """
    Defines an API exception
    """
    def __init__(self, message):
        super().__init__(message)

class ValidationException (Exception):
    """
    Defines an Validation exception
    """
    def __init__(self, message):
        super().__init__(message)


#def configure_logging(level):
#    """
#    Set logging level, based on the level names of the `logging` module.
#    - level (str): 'debug', 'info' or 'warning'
#    """
#    if level == 'debug':
#        level_name = logging.DEBUG
#    elif level == 'info':
#        level_name = logging.INFO
#    else:
#        level_name = logging.WARNING
#
#    logging.basicConfig(level=level_name,
#                        format='[%(asctime)s] %(levelname)s: %(message)s',
#                        datefmt='%Y-%m-%d@%H:%M:%S')


def parse_args():
    """
    Parses command-line options

    Raises:
        ArgumentException: if there are problems parsing arguments
    """
    parser = argparse.ArgumentParser(description="Utilities for QLIN", epilog="Ex: python qlin -u https://qlin-me-hybrid.staging.cqgc.hsj.rtss.qc.ca -t Termes_HPO_EXOG.txt.formatted.txt")
    parser.add_argument('-u', '--url', required=True, help="QLIN url, ex: https://qlin-me-hybrid.staging.cqgc.hsj.rtss.qc.ca")
    parser.add_argument('-m', '--mode', required=True, help="sequencting type either [germinal,somatic]")
    parser.add_argument('-t', '--termes', required=False, help="TermesHPO file for germinal mode, ex: Termes_HPO_EXOG.txt.formatted.txt")
    parser.add_argument('-a', '--analyses', required=False, help="Analyses file fot somatic mode, ex: Analyses.txt")
    parser.add_argument('-v', '--validate', action='store_true', required=False, help="Will validate input data with dry run")
    parser.add_argument('-f', '--fix', action='store_true', required=False, help="Will fix incorrect fields in Nanuq")
    parser.add_argument('-n', '--anonymise', action='store_true', required=False, help="Will anonymise patient infos")
    args = parser.parse_args()
    logger.info (args)
    if args.mode != 'germinal' and args.mode != 'somatic': raise ArgumentException("--mode must be one of [germinal,somatic]")
    if args.mode == 'germinal' and not args.termes: raise ArgumentException("Germinal mode does not have a termes argument")
    elif args.mode == 'somatic' and not args.analyses: raise ArgumentException("Somatic mode does not have an analyses argument")
    return(args)

def generate_random_string(length):
    """
    Utility function to generate a string of variable lenght using random ascii caracters

    Args:
        length (int)

    Returns:
        random_string (string): the resulting string

    Example:
        >>> generate_random_string(8)
        sherychh
    """
    characters = string.ascii_lowercase  # or string.ascii_letters, string.printable
    random_string = ''.join(random.choice(characters) for i in range(length))
    return random_string

class qlin:
    """
    Returns an object to interact with QLIN. The url and termes parameters are used to authenticate to the API.

    Args:
        url (string): the url to the QLIN API
    """

    def __init__(self, url):
        self.url                  = url
        self.authenticatedHeaders = self.authenticate(url)
#        configure_logging(level)

    def authenticate(self, url):
        """
        Authenticates a user to the QLIN API. Credentials will be fetched from ~/.qlin

        Args:
            url (string): the url to register

        Raises:
            APIException: if we fail to retreive a valid token

        Returns:
            authenticatedHeaders (object): the headers use to connect to the API service 
        """
    
        # loading variables from .env file
        qlin_config = dotenv_values(os.path.expanduser("~/.qlin"))
        authHeaders = {
            'Content-Type': 'application/x-www-form-urlencoded',
            'Accept': 'application/json'
        }
    
        # Enter your credentials here
        authPayload = {
            'email': qlin_config['email'],
            'password': qlin_config['password']
        }
        # Validate email and password
        if not authPayload['email'] or not authPayload['password']:
            raise APIException("Error: Email and password must be provided.")
    
        # Authentification
        response = requests.post(url+'/api/v1/auth/login', headers=authHeaders, data=authPayload)
    
        # Check if the request was successful
        if response.status_code == 200:
            # Extract the token
            token = response.json().get('token')
#            logging.debug(f"Authentication successful")
            authenticatedHeaders = {
              'Content-Type': 'application/json',
              'Accept': 'application/json',
              'Authorization': f'Bearer {token}'
            }
        else:
            raise APIException (f"Failed to retrieve token\nStatus code: {response.status_code}\n\nResponse:\n{json.dumps(response.json(),indent=2)}")

        return authenticatedHeaders
   

    def search_analysis(self, aliquot=None, sample=None, specimen=None, jhn=None, mrn=None, analysis_id=None, sequencing_id=None):
        """
        Implement method for endpoint '/api/v1/search/analysis'.
        Multiple arguments are treated in params as 'AND' condition.
        - aliquot: [str] CQGC lab id. Ex: '40250'
        # TODO: args for sample, specimen, mrn, jhn, ...
        - return : [list] Analysis, if found. 
                   Else, HTTP errors 400 (bad request) or 403 (forbidden).
        """
        endpoint = '/api/v1/search/analysis?'
        params   = ''
        if aliquot:     params += f'aliquot={aliquot}&'
        if analysis_id: params += f'analysis_id={analysis_id}&'
        if mrn:         params += f'mrn={mrn}&'

        response = requests.get(f"{self.url}{endpoint}{params}", headers=self.authenticatedHeaders)
        if response.status_code == 200:
            return response.json().get('analysis')
        else:
            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")

            
    def get_an_analysis(self, analysis_id):
        """
        GET analysis by `analysis_id` from Qlin
        - analysis_id: [str] ex: '822034'
        - return     : [dict] Analysis
        """
        endpoint = '/api/v1/analysis/'
        response = requests.get(f"{self.url}{endpoint}{analysis_id}", headers=self.authenticatedHeaders)
        if response.status_code == 200:
            return response.json()#.get('analysis')
        else:
            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\Analysis ID:\n{analysis_id}")

            
    def search_analysis_from_payload_all (self, analysis_payload):
        """
        Searches QLIN and return the analyses that matches all information from an analysis payload.

        Args:
            analysis_payload (json): the payload to search from

        Raises:
            a APIException if the response is not valid

        Returns:
            a list of analysis from QLIN
        """

        params = ""
        for patient in analysis_payload['patients']:
            if patient['jhn']: params += 'jhn=' + patient['jhn'] + '&'
            if patient['mrn']: params += 'mrn=' + patient['mrn'] + '&'
            params += 'aliquot=' + patient['aliquot'] + '&'
            params += 'sample=' + patient['sample'] + '&'
        params = params[:-1]
        response = requests.get(self.url + '/api/v1/search/analysis?' + str(params), headers=self.authenticatedHeaders)
        # Check if the request was successful
        if response.status_code == 200:
            response_json = response.json()
            return response_json
        else:
            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")


    def search_analysis_from_payload_jhn_mrn (self, analysis_payload):
        """
        Searches QLIN and return the analyses that matches the jhn and mrn included in an analysis payload.

        Args:
            analysis_payload (json): the payload to search from

        Raises:
            a APIException if the response is not valid

        Returns:
            a list of analysis from QLIN
        """

        params = ""
        for patient in analysis_payload['patients']:
            if patient['jhn']: params += 'jhn=' + patient['jhn'] + '&'
            if patient['mrn']: params += 'mrn=' + patient['mrn'] + '&'
        params = params[:-1]
        response = requests.get(self.url + '/api/v1/search/analysis?' + str(params), headers=self.authenticatedHeaders)
        # Check if the request was successful
        if response.status_code == 200:
            response_json = response.json()
            return response_json['analysis']
        else:
            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")



    def search_analysis_from_payload_aliquot_sample (self, analysis_payload):
        """
        Searches QLIN and return the analyses that matches the aliquot and sample id  included in an analysis payload.

        Args:
            analysis_payload (json): the payload to search from

        Raises:
            a APIException if the response is not valid

        Returns:
            a list of analysis from QLIN
        """

        params = ""
        for patient in analysis_payload['patients']:
            params += 'aliquot=' + patient['aliquot'] + '&'
            params += 'sample=' + patient['sample'] + '&'
        params = params[:-1]
        response = requests.get(self.url + '/api/v1/search/analysis?' + str(params), headers=self.authenticatedHeaders)
        # Check if the request was successful
        if response.status_code == 200:
            response_json = response.json()
            return response_json['analysis']
        else:
            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")


    def search_sequencing_from_aliquot_sample_in_payload (self, sequencing_payload):
        """
        Searches QLIN and return the sequncing that matches the aliquot and sample id included in an sequencing payload.

        Args:
            analysis_payload (json): the payload to search from

        Raises:
            a APIException if the response is not valid

        Returns:
            a list of analysis from QLIN
        """	
        params = ""
        for sequencing in sequencing_payload['sequencings']:
            params += 'aliquot=' + sequencing['aliquot'] + '&'
            params += 'sample=' + sequencing['sample'] + '&'
        params = params[:-1]
        response = requests.get(self.url + '/api/v1/search/analysis?' + str(params), headers=self.authenticatedHeaders)
        # Check if the request was successful
        if response.status_code == 200:
            response_json = response.json()
            return response_json['analysis']
        else:
            raise APIException (f"Failed search sequencing\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")


    def get_analyses_payloads_from_EXOG_Termes(self, file_termes):
        """
        Parses the content of the EXOG termes file into a dataframe, a manipulates the data structure into a valid analysis payload

        Args:
            file_termes (string): the name of the file to convert into a payload

        Returns:
            a list of analyses payloads (one per case)
        """
        # Parse using Panda read_csv utility
        termesDF_column_names = [ 'run_name', 'aliquot', 'Statut séquençage', 'sample', 'specimen', 'mrn', 'jhn', 'birth_date', 'sex',
                   'Design de famille', 'Famille ID',  'family_member', 'Statut', 'HPO', 'In current batch', 'Run associée',
                   'FAM', 'Individual', 'Father', 'Mother', 'Sex', 'Pheno', 'sequencing_types', 'type', 'analysis_code',
                   'priority', 'first_name', 'last_name', 'organization_id', 'laboratory_id' ]

        termesDF = pd.read_csv(file_termes, sep='\t', header=None, names=termesDF_column_names)
        
        # Replace known values
        mappings = {'Mâle': 'MALE', 'Femelle': 'FEMALE', 'Inconnu': 'UNKNOWN', 'Proband': 'PROBAND', 'Père': 'FTH', 'Mère':'MTH', 'Frère':'BRO', 'Soeur':'SIS' }
        termesDF = termesDF.replace(mappings)

        # Redefine column types after overwritten types from pd.read_csv
        for col in termesDF_column_names:
            if col not in termesDF.columns:
                termesDF[col] = pd.Series(dtype='str')
            else:
                termesDF[col] = termesDF[col].astype('str')

        # Iterate through cases to create analysis payload
        analyses_payloads = []
        for familyIDs in termesDF['FAM'].unique():
            analysis_payload = {}
            termesFamilyDF = termesDF[termesDF['FAM'] == familyIDs]
            termesFamilyProbandDF = termesFamilyDF[termesFamilyDF['family_member']=='PROBAND'].iloc[0]

            # Setting up defaulkt values for EXOG
            analysis_payload['type'] = "GERMLINE"
            analysis_payload['sequencing_types'] = ["WXS"]
            analysis_payload['diagnosis_hypothesis'] = 'unknown_hypothesis'

            # Update patients information using the Termes Dataframe
            patients=[]
            for index, individual in termesFamilyDF.iterrows():
                patient = {}
                patient['family_member'] = individual['family_member']
                if patient['family_member'] == 'MTH': patient['family_member'] = 'MOTHER'
                if patient['family_member'] == 'FTH': patient['family_member'] = 'FATHER'
                if (str(individual['jhn']) != 'nan'): patient['jhn'] = individual['jhn']
                patient['mrn'] = individual['mrn']
                patient['sex'] = individual['sex']
                patient['birth_date'] = individual['birth_date']
                patient['sample'] = individual['sample']
                patient['specimen'] = individual['specimen']
                patient['run_name'] = individual['run_name']                
                patient['aliquot'] = individual['aliquot']
                if patient['family_member'] != 'PROBAND':
                    patient['affected'] = True if individual['Statut'] == 'Affecté' else False
                    patient['status'] = 'NOW'
                if patient['family_member'] == 'PROBAND' or patient['affected'] == True:
                    if str(individual['HPO']) == 'nan': individual['HPO'] = 'HP:0000005'
                    clinical_signs = []
                    for hpo in str(individual['HPO']).split():
                        clinical_signs.append({'code': hpo.rstrip(';.,'), 'observed': True})
                    patient['clinical'] = {'signs':clinical_signs}
                patients.append(patient)
            analysis_payload['patients'] = patients
            analyses_payloads.append(analysis_payload)
        return analyses_payloads


    def get_analyses_payloads_from_EXOS_analyses(self, file_analyses):
        """
        Parses the content of the EXOS analyses.txt file into a dataframe, a manipulates the data structure into valid analysis payloads.

        Args:
            file_analyses (string): the name of the file to convert into a payload

        Raises:
            a ValidationException if the file_analyses is malformed

        Returns:
            a list of analyses payloads (one per case)
        """
        # Parse using Panda read_csv utility
        analysesDF_column_names = ['Date', 'Analysis_Run', 'Analysis_Sample', 'Normal_Run', 'Normal_Sample', 'Lab_Individual', 'Lab_Sample', 'Analysis_Name', 'Analysis_Type']
        analysesDF = pd.read_csv(file_analyses, sep='\t', header=None, names=analysesDF_column_names)

        # Iterate through cases to create analysis payload
        analyses_payloads = []
        for index, analyse in analysesDF.iterrows():
            patient = {}
            if analyse['Analysis_Type'] != 'WES_somatic-tumor_only':
                raise ValidationException(f"The analyses file for EXOS contains Analysis_Type different then WES_somatic-tumor_only")
            patient['family_member'] = 'PROBAND'
            patient['run_name'] = str(analyse['Analysis_Run'])
            patient['aliquot'] = str(analyse['Analysis_Sample'])
            patient['clinical'] = {'signs':[{'code':'HP:0000005', 'observed': True}]}
            analysis_payload = {
                'type':'SOMATIC_TUMOR_ONLY', 
                'sequencing_types':['WXS'], 
                'diagnosis_hypothesis':'unknown_hypothesis',
                'patients':[patient]
            }
            analyses_payloads.append(analysis_payload)
        return analyses_payloads


    def get_analyses_payloads_from_TRAN_analyses(self, file_analyses):
        """
        Parses the content of the TRAN analyses.txt file into a dataframe, a manipulates the data structure into valid analysis payloads

        Args:
            file_analyses (string): the name of the file to convert into a payload

        Raises:
            a ValidationException if the file_analyses is malformed

        Returns:
            a list of analyses payloads (one per case)
        """

        # Parse using Panda read_csv utility
        analysesDF_column_names = ['Date', 'Analysis_Run', 'Analysis_Sample', 'Normal_Run', 'Normal_Sample', 'Lab_Individual', 'Lab_Sample', 'Analysis_Name', 'Analysis_Type']
        analysesDF = pd.read_csv(file_analyses, sep='\t', header=None, names=analysesDF_column_names)

        # Iterate through cases to create analysis payload
        analyses_payloads = []
        for index, analyse in analysesDF.iterrows():
            patient = {}
            if analyse['Analysis_Type'] != 'RNASeq_somatic':
                raise ValidationException(f"The analyses file for TRAN contains Analysis_Type different then RNASeq_somatic")
            patient['family_member'] = 'PROBAND'
            patient['run_name'] = str(analyse['Analysis_Run'])
            patient['aliquot'] = str(analyse['Analysis_Sample'])
            patient['clinical'] = {'signs':[{'code':'HP:0000005', 'observed': True}]}
            analysis_payload = {
                'type':'SOMATIC_TUMOR_ONLY',
                'sequencing_types':['WTS'],
                'diagnosis_hypothesis':'unknown_hypothesis',
                'patients':[patient]
            }
            analyses_payloads.append(analysis_payload)
        return analyses_payloads


    def update_validate_EXOG_analysis_payload_from_nanuq(self, analysis_payload):
        """
        Compares and combine one EXOG case data from an analysis payload and from a Nanuq json object. Fixes common differences and fails raising a ValidationException otherwise.

        Args:
            analysis_payload (string): the name of the file to convert into a payload

        Raises:
            a ValidationException if miss matches exists between the payload and Nanuq

        Returns:
            an updated analysis payload to feed QLIN
        """
        nq = Nanuq()
        # Iterate trought family members to merge data from Nanuq and termes HPO
        for patient in analysis_payload['patients']:
            sample_nanuq = json.loads(nq.get_clinical_sample(patient['aliquot']))[0]
            if (str(sample_nanuq['labAliquotId']) != str(patient['aliquot'])):
                raise ValidationException(f"aliquot mismatch: {patient['aliquot']}, {sample_nanuq['labAliquotId']}")
            if (str(sample_nanuq['ldmSampleId']) != str(patient['sample']) ):
                raise ValidationException (f"ldmSampleId mismatch: {sample_nanuq['ldmSampleId']}, {patient['sample']}" )
            if (str(sample_nanuq['ldmSpecimenId']) != str(patient['specimen']) ):
                raise ValidationException (f"ldmSpecimenId mismatch: {sample_nanuq['ldmSpecimenId']}, {patient['specimen']}" )
            if (str(sample_nanuq['patient']['mrn']) != str(patient['mrn']) ):
                if (str(sample_nanuq['patient']['mrn']) != f"{int(patient['mrn']):08d}"):
                    raise ValidationException (f"mrn mismatch: cannot fix {sample_nanuq['patient']['mrn']}, {patient['mrn']}" )
                else:
                    logger.info (f"mrn mismatch: fixing {patient['mrn']} to {sample_nanuq['patient']['mrn']}" )
                    patient['mrn'] = '0'+str(patient['mrn'])
            try:
                if (str(sample_nanuq['patient']['ramq']) != str(patient['jhn']) ):
                    raise ValidationException (f"ramq mismatch: {sample_nanuq['patient']['ramq']}, {patient['jhn']}" )
            except KeyError:
                logging.info (f"Sample {patient['aliquot']} doesn't have ramq but mrn is valid")
            if (str(sample_nanuq['patient']['sex']) != str(patient['sex']) ):
                raise ValidationException (f"sex mismatch: {sample_nanuq['patient']['sex']}, {patient['sex']}" )

            #Fixing fields for MOTHER and FATHER
            if (str(sample_nanuq['patient']['familyMember']) == "MTH"): sample_nanuq['patient']['familyMember'] = 'MOTHER'
            if (str(sample_nanuq['patient']['familyMember']) == "FTH"): sample_nanuq['patient']['familyMember'] = 'FATHER'
            if (str(sample_nanuq['patient']['familyMember']) != str(patient['family_member']) ):
                raise ValidationException (f"familyMember mismatch: {sample_nanuq['patient']['familyMember']}, {patient['family_member']}" )

            # Fix improper panel codes from Nanuq
            if sample_nanuq['panelCode'] == 'PGDI':
                sample_nanuq['panelCode'] = 'EIDI'

           # Define values from nanuq
            patient['first_name'] = sample_nanuq['patient']['firstName']
            patient['last_name'] = sample_nanuq['patient']['lastName']
            patient['laboratory_id'] = sample_nanuq['ldm']
            patient['organization_id'] = sample_nanuq['patient']['ep']
            patient['sample_code'] = sample_nanuq['sampleType']
            patient['specimen_code'] = sample_nanuq['specimenType']
            if ( patient['family_member'] == 'PROBAND'):
                analysis_payload['analysis_code'] = sample_nanuq['panelCode']
                analysis_payload['priority'] = sample_nanuq['priority'] if (sample_nanuq['priority'] != '') else 'ROUTINE'

        return analysis_payload


    def update_validate_EXOS_TRAN_analysis_payload_from_nanuq(self, analysis_payload):
        """
        Compares and combine one EXOS or TRAN case data from an analysis payload and from a Nanuq json object.

        Args:
            analysis_payload (string): the name of the file to convert into a payload

        Returns:
            an updated analysis payload to feed QLIN
        """

        nq = Nanuq()
        # Iterate trought family members to merge data from Nanuq
        for patient in analysis_payload['patients']:
            sample_nanuq = json.loads(nq.get_clinical_sample(patient['aliquot']))[0]
            birth_date_nanuq = str(sample_nanuq['patient']['birthDate'])

            run = patient['run_name']
            m=re.search('(.*)/(.*)/(.*)', birth_date_nanuq)
            birth_date_clin = m.group(3) + "-" + m.group(2) + "-" + m.group(1)
            patient['birth_date'] = birth_date_clin
            patient['sample'] = str(sample_nanuq['ldmSampleId'])
            patient['specimen'] = str(sample_nanuq['ldmSpecimenId'])
            patient['mrn'] = str(sample_nanuq['patient']['mrn'])
            patient['jhn'] = str(sample_nanuq['patient']['ramq'])
            patient['sex'] = str(sample_nanuq['patient']['sex'])
            patient['first_name'] = str(sample_nanuq['patient']['firstName'])
            patient['last_name'] = str(sample_nanuq['patient']['lastName'])
            patient['laboratory_id'] = str(sample_nanuq['ldm'])
            patient['organization_id'] = str(sample_nanuq['patient']['ep'])
            patient['sample_code'] = str(sample_nanuq['sampleType'])
            patient['specimen_code'] = str(sample_nanuq['specimenType'])
            analysis_payload['analysis_code'] = str(sample_nanuq['panelCode'])
            analysis_payload['priority'] = str(sample_nanuq['priority']) if (sample_nanuq['priority'] != '') else 'ROUTINE'

        return analysis_payload


    def anonymise_analysis_payload(self, analysis_payload):
        """
        Returns a scrambled version of all sensible information in the analysis payload. This includes the mrn, jhn, first_name, last_name and sample.

        Args:
            analysis_payload (json): the payload to modify

        Returns:
            analysis_payload (json): the updated payload
        """
        for patient in analysis_payload['patients']:
             patient['mrn'] = generate_random_string(8)
             if 'jhn' in patient: patient['jhn'] = str(generate_random_string(4) + patient['jhn'][4:]).upper()
             patient['first_name'] = generate_random_string(8)
             patient['last_name'] = generate_random_string(8)
#             patient['sample'] = generate_random_string(8)
        return analysis_payload


    def validate_analysis_payload(self, analysis_payload):
        """
        Validates an analysis payload using QLIN '/api/v1/analysis' enpoint.

        Args:
            analysis_payload (json): json object containing the case analysis to create in QLIN

        Raises:
            ValidationException if a sample in the payload already exists in QLIN
            APIException if the API returns a status other then 200 (valid test)

        Returns:
            True if the analysis_payload validatese
        """
        analyses = self.search_analysis_from_payload_aliquot_sample(analysis_payload)
        if (len(analyses) > 0):
            raise ValidationException (f"Sample exists in QLIN:\n\nanalysis_payload: {analysis_payload}\n\nanalyses: {analyses}")

        response = requests.post(self.url + '/api/v1/analysis?validate-only=true', headers=self.authenticatedHeaders, data=json.dumps(analysis_payload))
        # Check if the request was successful
        if response.status_code != 200:
            raise APIException (f"Case creation could not validate\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(analysis_payload,indent=2)}")

        return True


    def push_analysis_payload(self, analysis_payload):
        """
        Creates patients and sequencings (the analysis) of a case by sending a POST to the '/api/v1/analysis' endpoint.

        Args:
            analysis_payload (json): json object containing the case analysis to create in QLIN

        Raises:
            APIException: if the API returns a status other 201 (insertion succes)

        Returns:
            analysis_payload (json): an updated payload containing the associated analysis_id, patient_id and sequencing_id created alongside the new case
        """

        response = requests.post(self.url + '/api/v1/analysis?validate-only=false', headers=self.authenticatedHeaders, data=json.dumps(analysis_payload))
        # Check if the request was successful
        if response.status_code == 201:
            response_json = response.json()
            for payload_patient in analysis_payload['patients']:
              for response_patient in response_json['patients']:
                if payload_patient['family_member'] == response_patient['family_member']:
                  payload_patient['analysis_id'] = response_json['analysis_id']
                  payload_patient['patient_id'] = response_patient['patient_id']
                  payload_patient['sequencing_id'] = response_patient['sequencings'][0]['sequencing_id']
            return analysis_payload
        else:
            raise APIException (f"Failed case creation\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(analysis_payload,indent=2)}")


    def get_sequencing_payload_EXOG(self, analysis_payload):
        """
        Generates the payload to pass to create_sequencings using the updated analysis_payload from create_analysis. Data consists of the sequencing info and the files generated by the bioinformatic analysis.

        Args:
            analysis_payload (json): a payload containing the case and the associated sequencing_ids

        Returns:
            sequencing_payload (json): a json payload
        """
        sequencing_payload = {}
        sequencings = []

        for patient in analysis_payload['patients']:
            run = patient['run_name']
            m=re.search('(.*)_(.*)_(.*)_(.*)', run)
            date_tmp=m.group(1)
            if len (date_tmp) == 8: #NovaX
                year = date_tmp[:4]
                if len(year) == 2: year = '20' + year
                date = year + "-" + date_tmp[4:6] + "-" + date_tmp[6:8]
            else: # Nova6000
                year = date_tmp[:2]
                if len(year) == 2: year = '20' + year
                date = year + "-" + date_tmp[2:4] + "-" + date_tmp[4:6]
            sequencer=m.group(2)
            run_id=m.group(3)
            flowcell=m.group(4)

            path_prefix = '/' + run + '_germinal/' + patient['aliquot']
            patient['resequencing'] = False
            files = [
                { "type": "ALIR",     "format": "CRAM", "path": path_prefix + ".dragen.WES_germinal.cram" },
                { "type": "ALIR",     "format": "CRAI", "path": path_prefix + ".dragen.WES_germinal.cram.crai" },
#                { "type": "ALIRPG"    "format": "CRAM", "path": path_prefix + ".masked.bam" },
#                { "type": "ALIRPG",   "format": "CRAI", "path": path_prefix + ".masked.bam.bai" },
#                { "type": "SNVPG",    "format": "VCF",  "path": path_prefix + ".maskedfiltered.norm.vcf.gz" },
#                { "type": "SNVPG",    "format": "TBI",  "path": path_prefix + ".maskedfiltered.norm.vcf.gz.tbi" },
#                { "type": "NRRV",     "format": "XLSX", "path": path_prefix + ".maskedfiltered.norm.VEP.annotated.xlsx" },
#                { "type": "PEC",      "format": "CSV",  "path": path_prefix + ".ploidy_estimation_metrics.csv" },
#                { "type": "ROH",      "format": "BED",  "path": path_prefix + ".roh.bed" },
                { "type": "SNV",      "format": "VCF",  "path": path_prefix + ".hard-filtered.gvcf.gz" },
                { "type": "SNV",      "format": "TBI",  "path": path_prefix + ".hard-filtered.gvcf.gz.tbi" },
                { "type": "GCNV",     "format": "VCF",  "path": path_prefix + ".cnv.vcf.gz" },
                { "type": "GCNV",     "format": "TBI",  "path": path_prefix + ".cnv.vcf.gz.tbi" },
                { "type": "GSV",      "format": "VCF",  "path": path_prefix + ".sv.vcf.gz" },
                { "type": "GSV",      "format": "TBI",  "path": path_prefix + ".sv.vcf.gz.tbi" },
                { "type": "SSUP",     "format": "TGZ",  "path": path_prefix + ".extra.tgz" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".seg.bw" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".hard-filtered.baf.bw" },
                { "type": "IGV",      "format": "BED",  "path": path_prefix + ".HyperExomeV2_combined_targets.bed" },
#                { "type": "IGV",      "format": "BED",  "path": path_prefix + ".KAPA_HyperExome_hg38_combined_targets.bed" },
#                { "type": "CNVVIS",   "format": "PNG",  "path": path_prefix + ".cnv.calls.png" },
                { "type": "COVGENE",  "format": "CSV",  "path": path_prefix + ".coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
                { "type": "QCRUN",    "format": "JSON", "path": path_prefix + ".QC_report.json" }
            ]
            if  patient['family_member'] == 'PROBAND':
                files.append( { "type": "NORM_VEP", "format": "VCF", "path": path_prefix + ".hard-filtered.formatted.norm.VEP.vcf.gz" } )
                files.append( { "type": "NORM_VEP", "format": "TBI", "path": path_prefix + ".hard-filtered.formatted.norm.VEP.vcf.gz.tbi" } )
#                if patient['organization_id'] == "CHUSJ":
                if len(patient['clinical']['signs']) > 0:
                    files.append( { "type": "EXOMISER", "format": "HTML", "path": path_prefix + ".exomiser.html" } )
                    files.append( { "type": "EXOMISER", "format": "JSON", "path": path_prefix + ".exomiser.json" } )
                    files.append( { "type": "EXOMISER", "format": "TSV",  "path": path_prefix + ".exomiser.variants.tsv" } )
            patient['files'] = files

            patient['experiment'] = {
                "platform": "Illumina",
                "sequencer": sequencer,
                "run_name": run,
                "run_date": date,
                "run_alias": run_id,
                "flowcell_id": flowcell,
                "is_paired_end": True,
                "fragment_size": 100,
                "experimental_strategy": "WXS",
                "capture_kit": "RocheKapaHyperExome",
#                "capture_kit": "RocheKapaHyperExomeV2",
                "bait_definition": "KAPA_HyperExome_hg38_capture_targets",
#                "bait_definition": "KAPA_HyperExomeV2_hg38_capture_targets",
                "protocol": "HyperPrep"
            }

            patient['workflow'] = { 
                "name": "Dragen", 
                "version": "4.4.4", 
                "genome_build": "GRCh38" 
            }

            sequencings.append(patient)

        sequencing_payload['sequencings'] = sequencings
        logger.debug(sequencing_payload)
        return sequencing_payload


    def get_sequencing_payload_EXOS(self, analysis_payload):
        """
        Generates the payload to pass to create_sequencings using the updated analysis_payload from create_analysis. Data consists of the sequencing info and the files generated by the bioinformatic analysis.

        Args:
            analysis_payload (json): a payload containing the case and the associated sequencing_ids

        Returns:
            sequencing_payload (json): a json payload
        """
        sequencing_payload = {}
        sequencings = []
        for patient in analysis_payload['patients']:
            run = patient['run_name']
            m=re.search('(.*)_(.*)_(.*)_(.*)', run)
            date_tmp=m.group(1)
            if len (date_tmp) == 8: #NovaX
                year = date_tmp[:4]
                if len(year) == 2: year = '20' + year
                date = year + "-" + date_tmp[4:6] + "-" + date_tmp[6:8]
            else: # Nova6000
                year = date_tmp[:2]
                if len(year) == 2: year = '20' + year
                date = year + "-" + date_tmp[2:4] + "-" + date_tmp[4:6]
            sequencer=m.group(2)
            run_id=m.group(3)
            flowcell=m.group(4)

            path_prefix = '/' + run + '_somatic/' + patient['aliquot']
            patient['resequencing'] = False
            files = [
                { "type": "ALIR",     "format": "CRAM", "path": path_prefix + ".dragen.WES_somatic-tumor_only.cram" },
                { "type": "ALIR",     "format": "CRAI", "path": path_prefix + ".dragen.WES_somatic-tumor_only.cram.crai" },
                { "type": "SNV",      "format": "VCF",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.hard-filtered.gvcf.gz" },
                { "type": "SNV",      "format": "TBI",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.hard-filtered.gvcf.gz.tbi" },
                { "type": "SCNV",     "format": "VCF",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.cnv.vcf.gz" },
                { "type": "SCNV",     "format": "TBI",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.cnv.vcf.gz.tbi" },
                { "type": "SSUP",     "format": "TGZ",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.extra.tgz" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".dragen.WES_somatic-tumor_only.seg.bw" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".dragen.WES_somatic-tumor_only.hard-filtered.baf.bw" },
                { "type": "IGV",      "format": "BED",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.KAPA_HyperExome_hg38_combined_targets.bed" },
                { "type": "COVGENE",  "format": "CSV",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
                { "type": "QCRUN",    "format": "JSON", "path": path_prefix + ".dragen.WES_somatic-tumor_only.QC_report.json" },
                { "type": "NORM_VEP", "format": "VCF",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.hard-filtered.norm.VEP.vcf.gz" },
                { "type": "NORM_VEP", "format": "TBI",  "path": path_prefix + ".dragen.WES_somatic-tumor_only.hard-filtered.norm.VEP.vcf.gz.tbi" }
            ]
            patient['files'] = files
            patient['experiment'] = {
                "platform": "Illumina",
                "sequencer": sequencer,
                "run_name": run,
                "run_date": date,
                "run_alias": run_id,
                "flowcell_id": flowcell,
                "is_paired_end": True,
                "fragment_size": 100,
                "experimental_strategy": "WXS",
                "capture_kit": "RocheKapaHyperExome",
                "bait_definition": "KAPA_HyperExome_hg38_capture_targets",
                "protocol": "HyperPlus"
            }
            patient['workflow'] = { "name": "Dragen", "version": "4.2.4", "genome_build": "GRCh38" }
            sequencings.append(patient)
        sequencing_payload['sequencings'] = sequencings
        return sequencing_payload


    def get_sequencing_payload_TRAN(self, analysis_payload):
        """
        Generates the payload to pass to create_sequencings using the updated analysis_payload from create_analysis. Data consists of the sequencing info and the files generated by the bioinformatic analysis.

        Args:
            analysis_payload (json): a payload containing the case and the associated sequencing_ids

        Returns:
            sequencing_payload (json): a json payload
        """
        sequencing_payload = {}
        sequencings = []
        for patient in analysis_payload['patients']:
            run = patient['run_name']
            m=re.search('(.*)_(.*)_(.*)_(.*)', run)
            date_tmp=m.group(1)
            if len (date_tmp) == 8: #NovaX
                year = date_tmp[:4]
                if len(year) == 2: year = '20' + year
                date = year + "-" + date_tmp[4:6] + "-" + date_tmp[6:8]
            else: # Nova6000
                year = date_tmp[:2]
                if len(year) == 2: year = '20' + year
                date = year + "-" + date_tmp[2:4] + "-" + date_tmp[4:6]
            sequencer=m.group(2)
            run_id=m.group(3)
            flowcell=m.group(4)

            path_prefix = '/' + run + '_somatic/' + patient['aliquot']
            patient['resequencing'] = False
            files = [
                { "type": "ALIR",       "format": "CRAM", "path": path_prefix + ".dragen.RNASeq_somatic_tumor.cram" },
                { "type": "ALIR",       "format": "CRAI", "path": path_prefix + ".dragen.RNASeq_somatic_tumor.cram.crai" },
                { "type": "SSNV",       "format": "VCF",  "path": path_prefix + ".dragen.RNASeq_somatic.hard-filtered.vcf.gz" },
                { "type": "SSNV",       "format": "TBI",  "path": path_prefix + ".dragen.RNASeq_somatic.hard-filtered.vcf.gz.tbi" },
                { "type": "SGENEXPSF",  "format": "TSV",  "path": path_prefix + ".dragen.RNASeq_somatic.quant.genes.sf" },
                { "type": "SOMFU",      "format": "VCF",  "path": path_prefix + ".dragen.RNASeq_somatic.fusion_candidates.vcf.gz" },
                { "type": "SOMFU",      "format": "TBI",  "path": path_prefix + ".dragen.RNASeq_somatic.fusion_candidates.vcf.gz.tbi" },
                { "type": "SRNASEQFCF", "format": "TSV",  "path": path_prefix + ".dragen.RNASeq_somatic.fusion_candidates.final" },
                { "type": "SRNASEQFCP", "format": "TSV",  "path": path_prefix + ".dragen.RNASeq_somatic.fusion_candidates.preliminary" },
                { "type": "COVGENE",    "format": "CSV",  "path": path_prefix + ".dragen.RNASeq_somatic.coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
                { "type": "EXP",        "format": "PNG",  "path": path_prefix + ".expression.png" },
                { "type": "FFUSION",    "format": "TSV",  "path": path_prefix + ".FusionSummary.tsv" },
                { "type": "FFUSION",    "format": "PDF",  "path": path_prefix + ".arriba.pdf" },
                { "type": "KMBCOR",     "format": "TXT",  "path": path_prefix + ".KM.BCOR-ITD-exon-15.fa.txt" },
                { "type": "KMFLT3",     "format": "TXT",  "path": path_prefix + ".KM.FLT3-ITD-exons-13-15.fa.txt" },
                { "type": "KMNPM1",     "format": "TXT",  "path": path_prefix + ".KM.NPM1-4ins-exons-10-11utr.fa.txt" },
                { "type": "KMUBTF",     "format": "TXT",  "path": path_prefix + ".KM.UBTF-ITD-2-exon-13.fa.txt" },
                { "type": "MDALL",      "format": "PDF",  "path": path_prefix + ".MD-ALL.umap.pdf" },
                { "type": "QCRUN",      "format": "JSON", "path": path_prefix + ".QC_report.json" },
                { "type": "SSUP",       "format": "TGZ",  "path": path_prefix + ".dragen.RNASeq_somatic.extra.tgz" },
            ]
            patient['files'] = files
            patient['experiment'] = {
                "platform": "Illumina",
                "sequencer": sequencer,
                "run_name": run,
                "run_date": date,
                "run_alias": run_id,
                "flowcell_id": flowcell,
                "is_paired_end": True,
                "fragment_size": 100,
                "experimental_strategy": "WTS",
                "capture_kit": "RocheKapaHyperExome",
                "bait_definition": "KAPA_HyperExome_hg38_capture_targets",
                "protocol": "HyperPlus"
            }
            patient['workflow'] = { "name": "Dragen", "version": "4.2.4", "genome_build": "GRCh38" }
            sequencings.append(patient)
        sequencing_payload['sequencings'] = sequencings
        return sequencing_payload


    def validate_sequencing_payload(self, sequencing_payload):
        """
        Validates a sequencing payload using QLIN '/api/v1/analysis/sequncings' enpoint.

        Args:
            sequencing_payload (json): json object containing the sequencing to associate to a case in QLIN

        Raises:
            APIException if the API returns a status other then 200 (valid test)

        Returns:
            True if the sequencing_payload validates
        """

        response = requests.post(self.url + '/api/v1/analysis/sequencings?validate-only=true', headers=self.authenticatedHeaders, data=json.dumps(sequencing_payload))
        # Check if the request was successful
        if response.status_code != 200:
            raise APIException (f"Sequencing creation could not validate\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(sequencing_payload,indent=2)}")

        return True


    def push_sequencing_payload(self, sequencing_payload):
        """
        Links the sequencings to the patients created previously in the electronic prescription by sending a POST to the '/api/v1/analysis/sequencings' endpoint.

        Args:
            sequencing_payload (json): json object containing the sequencings to associate in QLIN

        Raises:
            APIException: if the API returns a status other then 201 (insertion succes)

        Returns:
            response (json): the response of the API call
        """

        response = requests.post(self.url + '/api/v1/analysis/sequencings?validate-only=false', headers=self.authenticatedHeaders, data=json.dumps(sequencing_payload))
        # Check if the request was successful
#        if response.status_code == 201 or response.status_code == 500: 
        if response.status_code == 201:
### A enlever quand le bug response=500 sera corrige dans QLIN
#            print(f"Sequencing created {response.status_code}\n{json.dumps(sequencing_payload,indent=2)}\n{response.text}")
### A Remettre quand le bug response=500 sera corrige dans QLIN
#            return response.json()
            return response.json()
        else:
#            logging.error (f"Failed sequencing creation\n\n{response}")
            raise APIException (f"Failed sequencing creation\n\nStatus code: {response.status_code}\n\nResponse:\n{response}\n\nPayload:\n{json.dumps(sequencing_payload,indent=2)}")


    def get_pipeline_payload(self, sequencing_payload):
        """
        Generates the payload to pass to create_pipeline using the updated sequencing_payload from create_sequencings. Data consists of the sequencing ids to launch the pipeline.

        Args:
            sequencing_payload (json): a payload containing the sequencing ids to launch pipelines

        Returns:
            pipeline_payload (json): a json payload
        """

        pipeline_payload = {}
        sequencing_ids = []
        for sequencing in sequencing_payload['sequencings']:
            sequencing_ids.append(sequencing['sequencing_id'])
        pipeline_payload['sequencing_ids'] = sequencing_ids
        return pipeline_payload


    def validate_pipeline_payload(self, pipeline_payload):
        """
        Validates a pipeline payload using QLIN '/api/v1/analysis/run' enpoint.

        Args:
            pipeline_payload (json): json object containing the sequencing to associate to a case in QLIN

        Raises:
            APIException if the API returns a status other then 200 (valid test)

        Returns:
            True if the sequencing_payload validates
        """

        response = requests.post(self.url + '/api/v1/analysis/run?validate-only=true', headers=self.authenticatedHeaders, data=json.dumps(pipeline_payload))
        # Check if the request was successful
        if response.status_code != 200:
            raise ValidationException (f"Pipeline creation could not validate\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(pipeline_payload,indent=2)}")

        return True


    def push_pipeline_payload(self, pipeline_payload):
        """
        Triggers the pipeline that inserts sequencing data into QLIN by sending a POST to the '/api/v1/analysis/run' endpoint.

        Args:
            pipeline_payload (json): json object containing the sequencing ids to launch

        Raises:
            APIException: if the API returns a status other then 201 (insertion succes)

        Returns:
            response (json): the response of the API call
        """

        response = requests.post(self.url + '/api/v1/analysis/run?validate-only=false', headers=self.authenticatedHeaders, data=json.dumps(pipeline_payload))
        # Check if the request was successful
        if response.status_code == 201:
#            print(f"Pipeline started Case: {pipeline_payload}")
            return response.json()
        else:
            raise APIException (f"Failed pipeline start\n\nStatus code: {response.status_code}\n\nResponse:\n{response}\n\nPayload:\n{json.dumps(pipeline_payload,indent=2)}")


#    def load_analyses(self, file_analyses):
#        """
#        Parses the content of the analyses file into a valid dataframe
#
#        Args:
#            file_analyses (string): the name of the file to convert into dataframe
#
#        Returns:
#            analyses (pandas.DataFrame): a dataframe contaning the information in the analyses file
#        """
#        analysesDF_column_names = [
#            'Date',
#            'Analysis_Run',
#            'Analysis_Sample',
#            'Normal_Run',
#            'Normal_Sample',
#            'Lab_Individual',
#            'Lab_Sample',
#            'Analysis_Name',
#            'Analysis_Type'
#        ]
#        analysesDF = pd.read_csv(file_analyses, sep='\t', header=None, names=analysesDF_column_names)
#        return analysesDF
#
#    def reset_analyse_to_termes(self, analyse, run, termesDF):
#        """
#        Transforms the analyse dataframe and apprends it to the termeDF local object
#
#        Args:
#            analyse (json): Dataframe coming from the somatic exomes analyses.txt file
#
#        Returns:
#            The modified termesDF dataframe
#        """
#
#
#        nq = Nanuq()
#        sample_nanuq = json.loads(nq.get_sample(analyse['Analysis_Sample']))[0]
##        date_string_dd_mm_yyyy = "14/07/2025"
##        date_object = datetime.strptime(date_string_dd_mm_yyyy, "%d/%m/%Y")
##        date_string_yyyy_mm_dd = date_object.strftime("%Y-%m-%d")
#
#        formatted_birthdate = datetime.strptime(sample_nanuq['patient']['birthDate'], "%d/%m/%Y").strftime("%Y-%m-%d")
#        new_row_data = {
#            'run_name': run, 
#            'aliquot': sample_nanuq['labAliquotId'],
#            'sample': sample_nanuq['ldmSampleId'],
#            'specimen': sample_nanuq['ldmSpecimenId'],
#            'mrn': sample_nanuq['patient']['mrn'],
#            'birth_date': formatted_birthdate,
#            'sex': sample_nanuq['patient']['sex'],
#            'family_member': sample_nanuq['patient']['familyMember'],
#            'Statut': sample_nanuq['patient']['status'],
#            'sequencing_types': None,
#            'type': None,
#            'analysis_code': None,
#            'priority': None,
#            'first_name': None,
#            'last_name': None,
#            'organization_id': None,
#            'laboratory_id': None
#        }
#        if 'ramq' in sample_nanuq['patient']:
#            new_row_data['jhn'] = sample_nanuq['patient']['ramq']
#
#        termesDF = pd.DataFrame(columns=self.termesDF_column_types.keys()).astype(self.termesDF_column_types)
#        termesDF.loc[len(termesDF)] = new_row_data
#        return termesDF
#
#    
#    def get_somatic_sequencings_payload(self, analysis_payload):
#        """
#        Generates the payload to pass to create_sequencings using the updated analysis_payload from create_analysis. Data consists of the sequencing info and the files generated by the bioinformatic analysis.
#
#        Args:
#            analysis_payload (json): a payload containing the case and the associated sequencing_ids
#
#        Returns:
#            sequencing_payload (json): a json payload
#        """
#
#        run = self.run
#        m=re.search('(.*)_(.*)_(.*)_(.*)', run)
#        date_tmp=m.group(1)
#        year = date_tmp[:2]
#        if len(year) == 2: year = '20' + year
#        date = year + "-" + date_tmp[2:4] + "-" + date_tmp[4:6]
#        sequencer=m.group(2)
#        run_id=m.group(3)
#        flowcell=m.group(4)
#
#        sequencing_payload = {}
#        sequencings = []
#
#        for patient in analysis_payload['patients']:
#            path_prefix = '/' + run + '_somatic/' + patient['aliquot'] + '.dragen.WES_somatic-tumor_only'
#            sequencing = {}
#            sequencing['sequencing_id'] = patient['sequencing_id']
#            sequencing['resequencing'] = False
#            sequencing['laboratory_id'] = patient['laboratory_id']
#            sequencing['sample'] = patient['sample']
#            sequencing['specimen'] = patient['specimen']
#            sequencing['specimen_code'] = "TUMOR"
#            sequencing['sample_code'] = "DNA"
#            sequencing['aliquot'] = patient['aliquot']
#            files = [
#                { "type": "ALIR",     "format": "CRAM", "path": path_prefix + ".cram" },
#                { "type": "ALIR",     "format": "CRAI", "path": path_prefix + ".cram.crai" },
#                { "type": "CNVVIS",   "format": "PNG",  "path": path_prefix + ".cnv.png" },
#                { "type": "SCNV",     "format": "VCF",  "path": path_prefix + ".cnv.vcf.gz" },
#                { "type": "SCNV",     "format": "TBI",  "path": path_prefix + ".cnv.vcf.gz.tbi" },
#                { "type": "COVGENE",  "format": "CSV",  "path": path_prefix + ".coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
#                { "type": "SSUP",     "format": "TGZ",  "path": path_prefix + ".extra.tgz" },
#                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".hard-filtered.baf.bw" },
#                { "type": "SSNV",     "format": "VCF",  "path": path_prefix + ".hard-filtered.gvcf.gz" },
#                { "type": "SSNV",     "format": "TBI",  "path": path_prefix + ".hard-filtered.gvcf.gz.tbi" },
#                { "type": "NORM_VEP", "format": "VCF",  "path": path_prefix + ".hard-filtered.norm.VEP.vcf.gz" },
#                { "type": "NORM_VEP", "format": "TBI",  "path": path_prefix + ".hard-filtered.norm.VEP.vcf.gz.tbi" },
#                { "type": "IGV",      "format": "BED",  "path": path_prefix + ".KAPA_HyperExome_hg38_combined_targets.bed" },
#                { "type": "QCRUN",    "format": "JSON", "path": path_prefix + ".QC_report.json" },
#                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".seg.bw" },
#            ]
#            sequencing['files'] = files
#            sequencing['experiment'] = {
#                "platform": "Illumina",
#                "sequencer": sequencer,
#                "run_name": run,
#                "run_date": date,
#                "run_alias": run_id,
#                "flowcell_id": flowcell,
#                "is_paired_end": True,
#                "fragment_size": 100,
#                "experimental_strategy": "WXS",
#                "capture_kit": "RocheKapaHyperExome",
#                "bait_definition": "KAPA_HyperExome_hg38_capture_targets",
#                "protocol": "HyperPlus"
#            }
#            sequencing['workflow'] = { "name": "Dragen", "version": "4.2.4", "genome_build": "GRCh38" }
#            sequencings.append(sequencing)
#        sequencing_payload['sequencings'] = sequencings
#        return sequencing_payload


def main():
    """
    This is the main function of the standalone app.

    Args:
        --run -r: The name of the run. Ex: 250620_A00516_0688_AHGYCYDSXF
        --termes -t: The name of the file containing informations about cases. Usually named Termes_HPO_EXOG.txt.formatted.txt
        --url -u: the url to the QLIN API. Usually https://qlin-me-hybrid.cqgc.hsj.rtss.qc.ca or https://qlin-me-hybrid.staging.cqgc.hsj.rtss.qc.ca
        --validate -v: Flag that indicates we only want to validate input data without creating any new patients and cases.
    """

    args = parse_args()
    nq = Nanuq()
#    q = qlin(args.run, args.url, args.termes, args.analyses, args.mode)
    q = qlin(args.url)
#    run = args.run
    termes = args.termes
    analyses = args.analyses
    if args.mode == 'germinal':

        analysis_payload = q.termes_to_analysis_payload(args.termes, args.anonymise, args.fix)
#        print (analysis_payload)
        analysis_payload = q.nanuq_validate_analysis_payload(analysis_payload)

#        termesDF = q.load_termes(termes)
#        if args.validate:
#        # Validates the input data by iterating trough families 
#            for family in termesDF['FAM'].unique():
#                print (f"Validating fam {family}")
#                family_members = termesDF[termesDF['FAM'] == family]
#        
#                # Iterate trought family members
#                for index, individual in family_members.iterrows():
#                    print (f"Validating patient {individual}")
#                    sample_nanuq = json.loads(nq.get_sample(individual['aliquot']))[0]
#                    family_members.loc[index] = q.fix_nanuq_termes (sample_nanuq, individual)
#        
#                # Create analysis 
#                test_payload = q.get_analysis_payload(family_members, args.fix)
#                if args.anonymise: test_payload = q.anonymize_analysis_payload(test_payload)
#                print (f"Test payload: {test_payload,}")
#                q.create_analysis(test_payload, True)
#            print (f"test_payload: {test_payload}")
#            print (f"Validation succesfull")
#    
#        else:
#            # Iterate trough families to create analysis and sequencings
#            for family in termesDF['FAM'].unique():
#                family_members = termesDF[termesDF['FAM'] == family]
#        
#                # Iterate trought family members to merge data from Nanuq and termes HPO 
#                for index, individual in family_members.iterrows():
#                    sample_nanuq = json.loads(nq.get_sample(individual['aliquot']))[0]
#                    family_members.loc[index] = q.fix_nanuq_termes (sample_nanuq, individual)
#        
#                # Create and push analysis
#                analysis_payload = q.get_analysis_payload(family_members, args.fix)
#                if args.anonymise: analysis_payload = q.anonymize_analysis_payload(analysis_payload)
#                print (f"PRE analysis_payload: {analysis_payload}")
#                analysis_payload = q.create_analysis(analysis_payload)
#                print (f"POST analysis_payload: {analysis_payload}")
#                # Push sequencings data
#                sequencing_payload = q.get_germinal_sequencings_payload(analysis_payload)
#                print (f"PRE sequencing_payload: {sequencing_payload}")
#                sequencing_response = q.create_sequencings(sequencing_payload)
#                print (f"POST sequencing_payload: {sequencing_payload}") 
#       
#                # Start pipeline
#                pipeline_payload = q.get_pipeline_payload(sequencing_payload)
#                print (f"PRE pipeline_payload: {pipeline_payload}")
#                pipeline_response = q.create_pipeline(pipeline_payload)
#                print (f"POST pipeline_payload: {pipeline_payload}")


    elif args.mode == 'somatic':
#        termesDF = pd.DataFrame(columns=q.termesDF_column_types.keys()).astype(q.termesDF_column_types)
        analyses_payloads = q.get_analyses_payloads_EXOS(analyses)
        print (f"Avant: {analyses_payloads}")
        for analysis_payload in analyses_payloads:
            analysis_payload = q.update_validate_EXOS_analysis_payload_from_nanuq(analysis_payload)
            print (f"Apres: {analysis_payload}")
#        analysesDF = q.load_analyses(analyses)
#
#        for index, analyse in analysesDF.iterrows():
#            print (analyse)
#            termesDF = q.reset_analyse_to_termes(analyse, termesDF)
#
#            # Validates the input analyses by iterating members
#            family_members = termesDF
#
#            # Iterate trought family members
#            for index, individual in family_members.iterrows():
#                sample_nanuq = json.loads(nq.get_sample(individual['aliquot']))[0]
#                family_members.loc[index] = q.fix_nanuq_termes (sample_nanuq, individual)
#
#            # Create and push analysis
#            analysis_payload = q.get_analysis_payload(family_members, args.fix)
#            if args.validate: 
#                analysis_payload = q.create_analysis(analysis_payload, True)
#                print (f"Validation succesfull")
#            else: 
#                analysis_payload = q.create_analysis(analysis_payload, False)
#                # Push sequencings data
#                sequencing_payload = q.get_somatic_sequencings_payload(analysis_payload)
#                sequencing_response = q.create_sequencings(sequencing_payload)
#                # Start pipeline
#                pipeline_payload = q.get_pipeline_payload(sequencing_payload)
#                pipeline_response = q.create_pipeline(pipeline_payload)

    print (f"Completed")

if __name__ == '__main__':
    main()

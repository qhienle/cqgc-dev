#!/usr/bin/env python3

"""
Utilities for the QLIN API that creates analyses, patients and sequencings in the system

Import as a package from other modules:

>>> python
from qlin import qlin
q = qlin(url, termes)

As a standalone app, get the standard usage:

>>> shell
clin.py --help

To validate a set of input data:

>>> shell
qlin.py \\
  --termes Termes_HPO_EXOG.txt.formatted.txt \\
  --url    https://qlin-me-hybrid.staging.cqgc.hsj.rtss.qc.ca \\
  --validate

To create analyses, patients and sequencings in QLIN

>>> shell
qlin.py \\
  --termes Termes_HPO_EXOG.txt.formatted.txt \\
  --url    https://qlin-me-hybrid.staging.cqgc.hsj.rtss.qc.ca

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
#To format dates properly
from datetime import datetime
#importing necessary functions from dotenv library
from dotenv import load_dotenv, dotenv_values

# Imports from current lib
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from lib.nanuq import Nanuq

__version__ = "0.2"

class ArgumentException (Exception):
    """
    Defines an ArgumentException
    """
    def __init__(self, message):
        super().__init__(message)

class ParsingException (Exception):
    """
    Defines a Parsingexception
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
    print (args)
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
        termes (string): a file containing cases to insert into QLIN
    """

    def __init__(self, url):
        self.url                  = url
        self.authenticatedHeaders = self.authenticate(url)


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
            print(f"Authentication successful")
            authenticatedHeaders = {
              'Content-Type': 'application/json',
              'Accept': 'application/json',
              'Authorization': f'Bearer {token}'
            }
        else:
            raise APIException (f"Failed to retrieve token\nStatus code: {response.status_code}\n\nResponse:\n{json.dumps(response.json(),indent=2)}")
    
        return authenticatedHeaders
   


    def get_analyses_payloads_EXOG(self, file_termes):
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
            analysis_payload['diagnosis_hypothesis'] = 'diagnosis_hypothesis'

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


    def update_validate_analysis_payload_from_nanuq(self, analysis_payload):
        nq = Nanuq()
        # Iterate trought family members to merge data from Nanuq and termes HPO 
        for patient in analysis_payload['patients']:
            sample_nanuq = json.loads(nq.get_clinical_sample(patient['aliquot']))[0]
            if (str(sample_nanuq['labAliquotId']) != str(patient['aliquot'])):
                raise ParsingException(f"aliquot mismatch: {patient['aliquot']}, {sample_nanuq['labAliquotId']}")
            if (str(sample_nanuq['ldmSampleId']) != str(patient['sample']) ):
                raise ParsingException (f"ldmSampleId mismatch: {sample_nanuq['ldmSampleId']}, {patient['sample']}" )
            if (str(sample_nanuq['ldmSpecimenId']) != str(patient['specimen']) ):
                raise ParsingException (f"ldmSpecimenId mismatch: {sample_nanuq['ldmSpecimenId']}, {patient['specimen']}" )
            if (str(sample_nanuq['patient']['mrn']) != str(patient['mrn']) ):
                if (str(sample_nanuq['patient']['mrn']) != f"{int(patient['mrn']):08d}"):
                    raise ParsingException (f"mrn mismatch: cannot fix {sample_nanuq['patient']['mrn']}, {patient['mrn']}" )
                else:
                    print (f"mrn mismatch: fixing {patient['mrn']} to {sample_nanuq['patient']['mrn']}" )
                    patient['mrn'] = '0'+str(patient['mrn'])
            try:
                if (str(sample_nanuq['patient']['ramq']) != str(patient['jhn']) ):
                    raise ParsingException (f"ramq mismatch: {sample_nanuq['patient']['ramq']}, {patient['jhn']}" )
            except KeyError:
                print (f"Sample {patient['aliquot']} doesn't have ramq but mrn is valid")
            if (str(sample_nanuq['patient']['sex']) != str(patient['sex']) ):
                raise ParsingException (f"sex mismatch: {sample_nanuq['patient']['sex']}, {patient['sex']}" )

            #Fixing fields for MOTHER and FATHER
            if (str(sample_nanuq['patient']['familyMember']) == "MTH"): sample_nanuq['patient']['familyMember'] = 'MOTHER'
            if (str(sample_nanuq['patient']['familyMember']) == "FTH"): sample_nanuq['patient']['familyMember'] = 'FATHER'
            if (str(sample_nanuq['patient']['familyMember']) != str(patient['family_member']) ):
                raise ParsingException (f"familyMember mismatch: {sample_nanuq['patient']['familyMember']}, {patient['family_member']}" )

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


    def anonymise_analysis_payload(self, analysis_payload):
        for patient in analysis_payload['patients']:
             patient['mrn'] = generate_random_string(8)
             if 'jhn' in patient: patient['jhn'] = str(generate_random_string(4) + patient['jhn'][4:]).upper()
             patient['first_name'] = generate_random_string(8)
             patient['last_name'] = generate_random_string(8)
#             patient['sample'] = generate_random_string(8)
        return analysis_payload



#    def search_analysis_all (self, analysis_payload):
#        params = ""
#        for patient in analysis_payload['patients']:
#            if patient['jhn']: params += 'jhn=' + patient['jhn'] + '&'
#            if patient['mrn']: params += 'mrn=' + patient['mrn'] + '&'
#            params += 'aliquot=' + patient['aliquot'] + '&'
#            params += 'sample=' + patient['sample'] + '&'
#        params = params[:-1]
#        response = requests.get(self.url + '/api/v1/search/analysis?' + str(params), headers=self.authenticatedHeaders)
#        # Check if the request was successful
#        if response.status_code == 200:
#            response_json = response.json()
#            return response_json
#        else:
#            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")
#             
#
#    def search_analysis_jhn_mrn (self, analysis_payload):
#        params = ""
#        for patient in analysis_payload['patients']:
#            if patient['jhn']: params += 'jhn=' + patient['jhn'] + '&'
#            if patient['mrn']: params += 'mrn=' + patient['mrn'] + '&'
#        params = params[:-1]
#        response = requests.get(self.url + '/api/v1/search/analysis?' + str(params), headers=self.authenticatedHeaders)
#        # Check if the request was successful
#        if response.status_code == 200:
#            response_json = response.json()
#            return response_json['analysis']
#        else:
#            raise APIException (f"Failed search analyses\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nparams:\n{params}")


    def search_analysis_aliquot_sample (self, analysis_payload):
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

    def validate_analysis_payload(self, analysis_payload):
        analyses = self.search_analysis_aliquot_sample(analysis_payload)
        if (len(analyses) > 0):
            raise ValidationException (f"Sample exists in QLIN:\n\nanalysis_payload: {analysis_payload}\n\nanalyses: {analyses}")

        response = requests.post(self.url + '/api/v1/analysis?validate-only=true', headers=self.authenticatedHeaders, data=json.dumps(analysis_payload))
        # Check if the request was successful
        if response.status_code != 200:
            raise ValidationException (f"Case creation could not validate\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(analysis_payload,indent=2)}")


    def push_analysis_payload(self, analysis_payload):
        """
        Creates patients and sequencings (the analysis) of a case by sending a POST to the '/api/v1/analysis' endpoint.

        Args:
            analysis_payload (json): json object containing the case analysis to create in QLIN

        Raises:
            APIException: if the API returns a status other then 200 (valid test) or 201 (insertion succes)

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
#                { "type": "CNVVIS",   "format": "PNG",  "path": path_prefix + ".cnv.calls.png" },
                { "type": "COVGENE",  "format": "CSV",  "path": path_prefix + ".coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
                { "type": "QCRUN",    "format": "JSON", "path": path_prefix + ".QC_report.json" }
            ]
            if  patient['family_member'] == 'PROBAND':
                files.append( { "type": "NORM_VEP", "format": "VCF", "path": path_prefix + ".hard-filtered.formatted.norm.VEP.vcf.gz" } )
                files.append( { "type": "NORM_VEP", "format": "TBI", "path": path_prefix + ".hard-filtered.formatted.norm.VEP.vcf.gz.tbi" } )
                if patient['organization_id'] == "CHUSJ":
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
                "bait_definition": "KAPA_HyperExome_hg38_capture_targets",
                "protocol": "HyperPrep"
            }
            patient['workflow'] = { "name": "Dragen", "version": "4.4.4", "genome_build": "GRCh38" }
            sequencings.append(patient)
        sequencing_payload['sequencings'] = sequencings
        return sequencing_payload


    def search_sequencing_aliquot_sample (self, sequencing_payload):
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


    def validate_sequencing_payload(self, sequencing_payload):
        analyses = self.search_sequencing_aliquot_sample(sequencing_payload)
        if (len(analyses) > 0):
            raise ValidationException (f"Sample exists in QLIN:\n\nsequencing_payload: {sequencing_payload}\n\nanalyses: {analyses}")

        response = requests.post(self.url + '/api/v1/analysis/sequencings?validate-only=true', headers=self.authenticatedHeaders, data=json.dumps(sequencing_payload))
        # Check if the request was successful
        if response.status_code != 200:
            raise ValidationException (f"Sequencing creation could not validate\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(sequencing_payload,indent=2)}")


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
        if response.status_code == 201:
            print(f"Sequencing created {response.status_code}\n{json.dumps(sequencing_payload,indent=2)}\n{response.text}")
            return response.json()
        else:
            print (f"Failed sequencing creation\n\n{response}")
            raise APIException (f"Failed sequencing creation\n\nStatus code: {response.status_code}\n\nResponse:\n{response}\n\nPayload:\n{json.dumps(sequencing_payload,indent=2)}")


#    def get_pipeline_payload(self, sequencing_payload):
#        """
#        Generates the payload to pass to create_pipeline using the updated sequencing_payload from create_sequencings. Data consists of the sequencing ids to launch the pipeline.
#
#        Args:
#            sequencing_payload (json): a payload containing the sequencing ids to launch pipelines
#
#        Returns:
#            pipeline_payload (json): a json payload
#        """
#
#        pipeline_payload = {}
#        sequencing_ids = []
#        for sequencing in sequencing_payload['sequencings']:
#            sequencing_ids.append(sequencing['sequencing_id'])
#        pipeline_payload['sequencing_ids'] = sequencing_ids
#        return pipeline_payload
#
#
#    def validate_pipeline_payload(self, pipeline_payload):
#
#
#    def push_pipeline_payload(self, pipeline_payload):
#        """
#        Triggers the pipeline that inserts sequencing data into QLIN by sending a POST to the '/api/v1/analysis/run' endpoint.
#
#        Args:
#            pipeline_payload (json): json object containing the sequencing ids to launch
#
#        Raises:
#            APIException: if the API returns a status other then 201 (insertion succes)
#
#        Returns:
#            response (json): the response of the API call
#        """
#
#        response = requests.post(self.url + '/api/v1/analysis/run?validate-only=false', headers=self.authenticatedHeaders, data=json.dumps(pipeline_payload))
#        # Check if the request was successful
#        if response.status_code == 201:
#            print(f"Pipeline started Case: {pipeline_payload}")
#            return response.json()
#        else:
#            raise APIException (f"Failed pipeline start\n\nStatus code: {response.status_code}\n\nResponse:\n{response}\n\nPayload:\n{json.dumps(pipeline_payload,indent=2)}")









    def germinal_termes_to_analysis_payload(self, file_termes, fix):
        """
        NEW!
        """
        nq = Nanuq()

        termesDF_column_names = [ 'run_name', 'aliquot', 'Statut séquençage', 'sample', 'specimen', 'mrn', 'jhn', 'birth_date', 'sex', 'Design de famille', 'Famille ID',  'family_member', 'Statut', 'HPO', 'In current batch', 'Run associée', 'FAM', 'Individual', 'Father', 'Mother', 'Sex', 'Pheno', 'sequencing_types', 'type', 'analysis_code', 'priority', 'first_name', 'last_name', 'organization_id', 'laboratory_id' ]

        termesDF = pd.read_csv(file_termes, sep='\t', header=None, names=termesDF_column_names)

        # Redefine column types after overwritten types from pd.read_csv
        for col in termesDF_column_names:
            if col not in termesDF.columns:
                termesDF[col] = pd.Series(dtype='str')
            else:
                termesDF[col] = termesDF[col].astype('str')

        # Replance known values
        mappings = {'Mâle': 'MALE', 'Femelle': 'FEMALE', 'Inconnu': 'UNKNOWN', 'Proband': 'PROBAND', 'Père': 'FTH', 'Mère':'MTH', 'Frère':'BRO', 'Soeur':'SIS' }
        termesDF = termesDF.replace(mappings)

        print (termesDF)
        print (termesDF.to_json())
        print (termesDF.to_json(orient='records', indent=4))
#        for family in termesDF['FAM'].unique():
#            family_members = termesDF[termesDF['FAM'] == family]
#    
#            # Iterate trought family members to merge data from Nanuq and termes HPO 
#            for index, individual in family_members.iterrows():
#                sample_nanuq = json.loads(nq.get_sample(individual['aliquot']))[0]
#                print (f"sample_nanuq {sample_nanuq}")
##                print (f"{nq.get_sample(41307)}")
#                family_members.loc[index] = self.fix_nanuq_termes (sample_nanuq, individual)
#    
#            # Create and push analysis
#            analysis_payload = self.get_analysis_payload(family_members, fix)
#            if anonymise: analysis_payload = self.anonymize_analysis_payload(analysis_payload)
        analysis_payload = self.nanuq_validate_analysis_payload(termesDF)
        #print (f"analysis_payload: {analysis_payload}")
#        if anonymise: analysis_payload = self.anonymize_analysis_payload(analysis_payload)
        return analysis_payload

    def load_termes(self, file_termes):
        """
        Parses the content of the termes file into a valid dataframe

        Args:
            file_termes (string): the name of the file to convert into dataframe

        Returns:
            termesDF (pandas.DataFrame): a dataframe contaning the information in the termes file
        """

        termesDF = pd.read_csv(file_termes, sep='\t', header=None, names=self.termesDF_column_types.keys(), dtype=self.termesDF_column_types)

        # Replance known values
        mappings = {'Mâle': 'MALE', 'Femelle': 'FEMALE', 'Inconnu': 'UNKNOWN', 'Proband': 'PROBAND', 'Père': 'FTH', 'Mère':'MTH', 'Frère':'BRO', 'Soeur':'SIS' }
        termesDF = termesDF.replace(mappings)

        # Create columns not present Termes for payload
#        termesDF['sequencing_types'] = None
#        termesDF['type'] = None
#        termesDF['analysis_code'] = None
#        termesDF['priority'] = None
#        termesDF['first_name'] = None
#        termesDF['last_name'] = None
#        termesDF['organization_id'] = None
#        termesDF['laboratory_id'] = None

        return termesDF

    def load_analyses(self, file_analyses):
        """
        Parses the content of the analyses file into a valid dataframe

        Args:
            file_analyses (string): the name of the file to convert into dataframe

        Returns:
            analyses (pandas.DataFrame): a dataframe contaning the information in the analyses file
        """
        analysesDF_column_names = [
            'Date',
            'Analysis_Run',
            'Analysis_Sample',
            'Normal_Run',
            'Normal_Sample',
            'Lab_Individual',
            'Lab_Sample',
            'Analysis_Name',
            'Analysis_Type'
        ]
        analysesDF = pd.read_csv(file_analyses, sep='\t', header=None, names=analysesDF_column_names)
        return analysesDF


    def nanuq_validate_analysis_payload(self, analysis_payload):
        nq = Nanuq()
        print (analysis_payload)
        for patient in analysis_payload['patients']:
            print (f"patient: {patient}")
            sample_nanuq = json.loads(nq.get_sample(patient['aliquot']))[0]
            print (f"sample_nanuq: {sample_nanuq}")

    def fix_nanuq_termes_old(self, nanuqJSON, termesDF):
        """
        Compares and combine one individual patient data from a Nanuq json object and the termes DataFrame. Fixes common differences and fails raising a ParsingException otherwise 

        Args:
            nanuqJSON (json): json object containing the individual clinical attributes submitted in Nanuq
            termesDF (pandas.DataFrame): the dataframe containing the patient info included in the termes file

        Raises:
            ParsingException

        Returns:
            termesDF (pandas.DataFrame): a modified dataframe contaning the aggregate of the json and the dataframe data
        """

        # JSON format
        #{
        #    "labAliquotId": "36272",
        #    "labAliquotReceptionDate": "16/05/2025",
        #    "labAliquotSubmissionDate": "16/05/2025",
        #    "ldm": "LDM-CHUSJ",
        #    "ldmSampleId": "DM253501",
        #    "ldmServiceRequestId": ".",
        #    "ldmSpecimenId": "ZR050933",
        #    "libType": "KAPA HyperExome",
        #    "panelCode": "TUHEM",
        #    "patient": {
        #        "birthDate": "18/04/2006",
        #        "ep": "CHUSJ",
        #        "familyMember": "PROBAND",
        #        "fetus": false,
        #        "firstName": "MICHAEL",
        #        "lastName": "METIVIER",
        #        "mrn": "02310192",
        #        "ramq": "METM06041815",
        #        "sex": "MALE",
        #        "status": "AFF"
        #    },
        #    "priority": "ROUTINE",
        #    "projectGroup": "Exome_Germinal",
        #    "projectName": "Exome_Germinal_CHUSJ",
        #    "sampleType": "DNA",
        #    "specimenType": "NBL"
        #}

        # TermesDF format
        #run_name             250530_A00516_0683_AH5FYMDMX2
        #aliquot                                      36272
        #Statut séquençage                         séquencé
        #sample                                    DM253501
        #specimen                                  ZR050933
        #mrn                                        2310192
        #ramq                                  METM06041815
        #birth_date                              2006-04-18
        #sex                                           Mâle
        #Design de famille                                -
        #Famille ID                                 2310192
        #family_member                              Proband
        #Statut                                     Affecté
        #HPO                                     HP:0006727
        #In current batch                                 X
        #Run associée                                   NaN
        #FAM                                        2310192
        #Individual                                   36272
        #Father                                           0
        #Mother                                           0
        #Sex                                              1
        #Pheno                                            2

        if (str(nanuqJSON['labAliquotId']) != str(termesDF['aliquot'])):
            raise ParsingException(f"aliquot mismatch: {termesDF['aliquot']}, {nanuqJSON['labAliquotId']}")
        if (str(nanuqJSON['ldmSampleId']) != str(termesDF['sample']) ):
            raise ParsingException (f"ldmSampleId mismatch: {nanuqJSON['ldmSampleId']}, {termesDF['sample']}" )
        if (str(nanuqJSON['ldmSpecimenId']) != str(termesDF['specimen']) ):
            raise ParsingException (f"ldmSpecimenId mismatch: {nanuqJSON['ldmSpecimenId']}, {termesDF['specimen']}" )
        if (str(nanuqJSON['patient']['mrn']) != str(termesDF['mrn']) ):
            if (str(nanuqJSON['patient']['mrn']) != f"{int(termesDF['mrn']):08d}"):
#            if (str(nanuqJSON['patient']['mrn']) != str('0'+str(termesDF['mrn'])) ):
                raise ParsingException (f"mrn mismatch: cannot fix {nanuqJSON['patient']['mrn']}, {termesDF['mrn']}" )
            else:
                print (f"mrn mismatch: fixing {termesDF['mrn']} to {nanuqJSON['patient']['mrn']}" )
                termesDF['mrn'] = '0'+termesDF['mrn']
        try:
            if (str(nanuqJSON['patient']['ramq']) != str(termesDF['jhn']) ):
                raise ParsingException (f"ramq mismatch: {nanuqJSON['patient']['ramq']}, {termesDF['jhn']}" )
        except KeyError:
            print (f"Sample {termesDF['aliquot']} doesn't have ramq but mrn is valid")
        if (str(nanuqJSON['patient']['sex']) != str(termesDF['sex']) ):
            raise ParsingException (f"sex mismatch: {nanuqJSON['patient']['sex']}, {termesDF['sex']}" )
        if (str(nanuqJSON['patient']['familyMember']) != str(termesDF['family_member']) ):
            raise ParsingException (f"familyMember mismatch: {nanuqJSON['patient']['familyMember']}, {termesDF['family_member']}" )

#        # Analysis/Sequencing Type
##        if nanuqJSON["projectGroup"] == "Exome_Somatique":
#        if nanuqJSON["libType"] == "KAPA HyperExome" and nanuqJSON["specimenType"] == "TUMOR":
#            termesDF['type'] = "SOMATIC_TUMOR_ONLY"
#            termesDF['sequencing_types'] = ["WXS"]
##        elif nanuqJSON["projectGroup"] == "Exome_Germinal":
#        elif nanuqJSON["libType"] == "KAPA HyperExome" and nanuqJSON["specimenType"] == "NBL":
#            termesDF['type'] = "GERMLINE"
#            termesDF['sequencing_types'] = ["WXS"]
#        elif nanuqJSON["libType"] == "Stranded Total RNA" and nanuqJSON["specimenType"] == "TUMOR":
#            termesDF['type'] = "SOMATIC_TUMOR_ONLY"
#            termesDF['sequencing_types'] = ["WTS"]
#        else:
#            raise ParsingException (f"Unknown sequence type {sample_nanuq} for {args.sample}")

        # Define values from nanuq

        # Fix improper panel codes from Nanuq
        if nanuqJSON['panelCode'] == 'PGDI':
            nanuqJSON['panelCode'] = 'EIDI'

        termesDF['analysis_code'] = nanuqJSON['panelCode']
        termesDF['priority'] = nanuqJSON['priority'] if (nanuqJSON['priority'] != '') else 'ROUTINE'
        termesDF['first_name'] = nanuqJSON['patient']['firstName']
        termesDF['last_name'] = nanuqJSON['patient']['lastName']
        termesDF['laboratory_id'] = nanuqJSON['ldm']
        termesDF['organization_id'] = nanuqJSON['patient']['ep']

        return termesDF

    def reset_analyse_to_termes(self, analyse, run, termesDF):
        """
        Transforms the analyse dataframe and apprends it to the termeDF local object

        Args:
            analyse (json): Dataframe coming from the somatic exomes analyses.txt file

        Returns:
            The modified termesDF dataframe
        """


        nq = Nanuq()
        sample_nanuq = json.loads(nq.get_sample(analyse['Analysis_Sample']))[0]
#        date_string_dd_mm_yyyy = "14/07/2025"
#        date_object = datetime.strptime(date_string_dd_mm_yyyy, "%d/%m/%Y")
#        date_string_yyyy_mm_dd = date_object.strftime("%Y-%m-%d")

        formatted_birthdate = datetime.strptime(sample_nanuq['patient']['birthDate'], "%d/%m/%Y").strftime("%Y-%m-%d")
        new_row_data = {
            'run_name': run, 
            'aliquot': sample_nanuq['labAliquotId'],
            'sample': sample_nanuq['ldmSampleId'],
            'specimen': sample_nanuq['ldmSpecimenId'],
            'mrn': sample_nanuq['patient']['mrn'],
            'birth_date': formatted_birthdate,
            'sex': sample_nanuq['patient']['sex'],
            'family_member': sample_nanuq['patient']['familyMember'],
            'Statut': sample_nanuq['patient']['status'],
            'sequencing_types': None,
            'type': None,
            'analysis_code': None,
            'priority': None,
            'first_name': None,
            'last_name': None,
            'organization_id': None,
            'laboratory_id': None
        }
        if 'ramq' in sample_nanuq['patient']:
            new_row_data['jhn'] = sample_nanuq['patient']['ramq']

        termesDF = pd.DataFrame(columns=self.termesDF_column_types.keys()).astype(self.termesDF_column_types)
        termesDF.loc[len(termesDF)] = new_row_data
        return termesDF

    def get_analysis_payload_from_termes(self, termesFAM, fix):
        """
        Generates the payload to pass to QLIN '/api/v1/analysis' endpoint using the case termes DataFrame

        Args:
            termesFAM (pandas.Dataframe): a dataframe containing the data for all members of the case
            fix (boolean): indicates if we modify patient names, accepting only alphabetic caracters

        Returns:
            analysis_payload (json): a json payload
        """

        analysis_payload = {}
        proband = termesFAM[termesFAM['family_member']=='PROBAND'].iloc[0]
        analysis_payload['type'] = proband['type']
        analysis_payload['analysis_code'] = proband['analysis_code']
#        if analysis_payload['analysis_code'] == 'VEOIB' : analysis_payload['analysis_code'] = 'VEOIBD'
        analysis_payload['priority'] = proband['priority']
        analysis_payload['diagnosis_hypothesis'] = 'diagnosis_hypothesis'
        analysis_payload['sequencing_types'] = proband['sequencing_types']
        patients=[]
        for index, individual in termesFAM.iterrows():
            patient = {}
            patient['family_member'] = individual['family_member']
            if patient['family_member'] == 'MTH': patient['family_member'] = 'MOTHER'
            if patient['family_member'] == 'FTH': patient['family_member'] = 'FATHER'
            patient['first_name'] = individual['first_name']
            regex = r'[^a-zA-Z]'
            if fix and re.findall(regex, patient['first_name']):
                first_name = patient['first_name']
                patient['first_name'] = re.sub(regex, '', patient['first_name'])
                print ("Patient {patient} first name {first_name} renamed to {patient['first_name']}")
            patient['last_name'] = individual['last_name']
            if fix and re.findall(regex, patient['last_name']):
                last_name = patient['last_name']
                patient['last_name'] = re.sub(regex, '', patient['last_name'])
                print ("Patient {patient} last name {last_name} renamed to {patient['last_name']}")
            if (str(individual['jhn']) != 'nan'): patient['jhn'] = individual['jhn']
            patient['mrn'] = individual['mrn']
            patient['sex'] = individual['sex']
            patient['birth_date'] = individual['birth_date']
            patient['organization_id'] = individual['organization_id']
            patient['laboratory_id'] = individual['laboratory_id']
            patient['sample'] = individual['sample']
            patient['specimen'] = individual['specimen']
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
        return analysis_payload

    def anonymize_analysis_payload (self, analysis_payload):
        """
        Returns a scrambled version of all sensible information in the analysis payload. This includes the mrn, jhn, first_name and last_name.

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
        return analysis_payload


    def create_analysis (self, analysis_payload, validatation=False):
        """
        Creates patients and sequencings (the analysis) of a case by sending a POST to the '/api/v1/analysis' endpoint.

        Args:
            analysis_payload (json): json object containing the case analysis to create in QLIN
            validatation (Boolean): a flag that indicates this is a dry run (no insertion in QLIN).

        Raises:
            APIException: if the API returns a status other then 200 (valid test) or 201 (insertion succes)

        Returns:
            analysis_payload (json): an updated payload containing the associated analysis_id, patient_id and sequencing_id created alongside the new case
        """

        response = requests.post(self.url + '/api/v1/analysis?validate-only=' + str(validatation), headers=self.authenticatedHeaders, data=json.dumps(analysis_payload))
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
    
        elif response.status_code == 200:
            return None
    
        else:
            raise APIException (f"Failed case creation\n\nStatus code: {response.status_code}\n\nResponse:\n{response.text}\n\nPayload:\n{json.dumps(analysis_payload,indent=2)}")
    
    
    def get_germinal_sequencings_payload(self, analysis_payload):
        """
        Generates the payload to pass to create_sequencings using the updated analysis_payload from create_analysis. Data consists of the sequencing info and the files generated by the bioinformatic analysis.

        Args:
            analysis_payload (json): a payload containing the case and the associated sequencing_ids

        Returns:
            sequencing_payload (json): a json payload
        """

        run = self.run
        m=re.search('(.*)_(.*)_(.*)_(.*)', run)
        date_tmp=m.group(1)
        year = date_tmp[:2]
        if len(year) == 2: year = '20' + year
        date = year + "-" + date_tmp[2:4] + "-" + date_tmp[4:6]
        sequencer=m.group(2)
        run_id=m.group(3)
        flowcell=m.group(4)
    
        sequencing_payload = {}
        sequencings = []
        for patient in analysis_payload['patients']:
            path_prefix = '/' + run + '_germinal/' + patient['aliquot']
            sequencing = {}
            sequencing['sequencing_id'] = patient['sequencing_id']
            sequencing['resequencing'] = False
            sequencing['laboratory_id'] = patient['laboratory_id']
            sequencing['sample'] = patient['sample']
            sequencing['specimen'] = patient['specimen']
            sequencing['specimen_code'] = "NBL"
            sequencing['sample_code'] = "DNA"
            sequencing['aliquot'] = patient['aliquot']
            files = [
                { "type": "ALIR",     "format": "CRAM", "path": path_prefix + ".cram" },
                { "type": "ALIR",     "format": "CRAI", "path": path_prefix + ".cram.crai" },
                { "type": "SNV",      "format": "VCF",  "path": path_prefix + ".hard-filtered.gvcf.gz" },
                { "type": "SNV",      "format": "TBI",  "path": path_prefix + ".hard-filtered.gvcf.gz.tbi" },
                { "type": "GCNV",     "format": "VCF",  "path": path_prefix + ".cnv.vcf.gz" },
                { "type": "GCNV",     "format": "TBI",  "path": path_prefix + ".cnv.vcf.gz.tbi" },
                { "type": "GSV",      "format": "VCF",  "path": path_prefix + ".sv.vcf.gz" },
                { "type": "GSV",      "format": "TBI",  "path": path_prefix + ".sv.vcf.gz.tbi" },
                { "type": "SSUP",     "format": "TGZ",  "path": path_prefix + ".extra.tgz" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".seg.bw" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".hard-filtered.baf.bw" },
                { "type": "IGV",      "format": "BED",  "path": path_prefix + ".KAPA_HyperExome_hg38_combined_targets.bed" },
                { "type": "CNVVIS",   "format": "PNG",  "path": path_prefix + ".cnv.calls.png" },
                { "type": "COVGENE",  "format": "CSV",  "path": path_prefix + ".coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
                { "type": "QCRUN",    "format": "JSON", "path": path_prefix + ".QC_report.json" }
            ]
            if  patient['family_member'] == 'PROBAND':
                files.append( { "type": "NORM_VEP", "format": "VCF", "path": path_prefix + ".hard-filtered.formatted.norm.VEP.vcf.gz" } )
                files.append( { "type": "NORM_VEP", "format": "TBI", "path": path_prefix + ".hard-filtered.formatted.norm.VEP.vcf.gz.tbi" } )
                if patient['organization_id'] == "CHUSJ":
                    files.append( { "type": "EXOMISER", "format": "HTML", "path": path_prefix + ".exomiser.html" } )
                    files.append( { "type": "EXOMISER", "format": "JSON", "path": path_prefix + ".exomiser.json" } )
                    files.append( { "type": "EXOMISER", "format": "TSV",  "path": path_prefix + ".exomiser.variants.tsv" } )
            sequencing['files'] = files 
            sequencing['experiment'] = {
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
                "protocol": "HyperPrep"
            }
            sequencing['workflow'] = { "name": "Dragen", "version": "4.2.4", "genome_build": "GRCh38" }
            sequencings.append(sequencing)
        sequencing_payload['sequencings'] = sequencings
        return sequencing_payload

   
    def get_somatic_sequencings_payload(self, analysis_payload):
        """
        Generates the payload to pass to create_sequencings using the updated analysis_payload from create_analysis. Data consists of the sequencing info and the files generated by the bioinformatic analysis.

        Args:
            analysis_payload (json): a payload containing the case and the associated sequencing_ids

        Returns:
            sequencing_payload (json): a json payload
        """

        run = self.run
        m=re.search('(.*)_(.*)_(.*)_(.*)', run)
        date_tmp=m.group(1)
        year = date_tmp[:2]
        if len(year) == 2: year = '20' + year
        date = year + "-" + date_tmp[2:4] + "-" + date_tmp[4:6]
        sequencer=m.group(2)
        run_id=m.group(3)
        flowcell=m.group(4)

        sequencing_payload = {}
        sequencings = []

        for patient in analysis_payload['patients']:
            path_prefix = '/' + run + '_somatic/' + patient['aliquot'] + '.dragen.WES_somatic-tumor_only'
            sequencing = {}
            sequencing['sequencing_id'] = patient['sequencing_id']
            sequencing['resequencing'] = False
            sequencing['laboratory_id'] = patient['laboratory_id']
            sequencing['sample'] = patient['sample']
            sequencing['specimen'] = patient['specimen']
            sequencing['specimen_code'] = "TUMOR"
            sequencing['sample_code'] = "DNA"
            sequencing['aliquot'] = patient['aliquot']
            files = [
                { "type": "ALIR",     "format": "CRAM", "path": path_prefix + ".cram" },
                { "type": "ALIR",     "format": "CRAI", "path": path_prefix + ".cram.crai" },
                { "type": "CNVVIS",   "format": "PNG",  "path": path_prefix + ".cnv.png" },
                { "type": "SCNV",     "format": "VCF",  "path": path_prefix + ".cnv.vcf.gz" },
                { "type": "SCNV",     "format": "TBI",  "path": path_prefix + ".cnv.vcf.gz.tbi" },
                { "type": "COVGENE",  "format": "CSV",  "path": path_prefix + ".coverage_by_gene.GENCODE_CODING_CANONICAL.csv" },
                { "type": "SSUP",     "format": "TGZ",  "path": path_prefix + ".extra.tgz" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".hard-filtered.baf.bw" },
                { "type": "SSNV",     "format": "VCF",  "path": path_prefix + ".hard-filtered.gvcf.gz" },
                { "type": "SSNV",     "format": "TBI",  "path": path_prefix + ".hard-filtered.gvcf.gz.tbi" },
                { "type": "NORM_VEP", "format": "VCF",  "path": path_prefix + ".hard-filtered.norm.VEP.vcf.gz" },
                { "type": "NORM_VEP", "format": "TBI",  "path": path_prefix + ".hard-filtered.norm.VEP.vcf.gz.tbi" },
                { "type": "IGV",      "format": "BED",  "path": path_prefix + ".KAPA_HyperExome_hg38_combined_targets.bed" },
                { "type": "QCRUN",    "format": "JSON", "path": path_prefix + ".QC_report.json" },
                { "type": "IGV",      "format": "BW",   "path": path_prefix + ".seg.bw" },
            ]
            sequencing['files'] = files
            sequencing['experiment'] = {
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
            sequencing['workflow'] = { "name": "Dragen", "version": "4.2.4", "genome_build": "GRCh38" }
            sequencings.append(sequencing)
        sequencing_payload['sequencings'] = sequencings
        return sequencing_payload


 
    def create_sequencings(self, sequencing_payload):
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
        if response.status_code == 201:
    #        print(f"Sequencing created {response.status_code}\n{json.dumps(sequencing_payload,indent=2)}\n{response.text}")
            return response.json()
        else:
            raise APIException (f"Failed sequencing creation\n\nStatus code: {response.status_code}\n\nResponse:\n{json.dumps(response.json(),indent=2)}\n\nPayload:\n{json.dumps(sequencing_payload,indent=2)}")

 
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
    
    
    def create_pipeline(self, pipeline_payload):
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
            print(f"Pipeline started Case: {pipeline_payload}")
            return response.json()
        else:
            raise APIException (f"Failed pipeline start\n\nStatus code: {response.status_code}\n\nResponse:\n{response}\n\nPayload:\n{json.dumps(pipeline_payload,indent=2)}")

def main():
    """
    This is the main function of the standalone app.

    Args:
        --run -r: The name of the run. Ex: 250620_A00516_0688_AHGYCYDSXF
        --termes -t: The name of the file containing informations about cases. Usually named Termes_HPO_EXOG.txt.formatted.txt
        --url -u: the url to the QLIN API. Usually https://qlin-me-hybrid.cqgc.hsj.rtss.qc.ca or https://qlin-me-hybrid.staging.cqgc.hsj.rtss.qc.ca
        --validate -v: Flag that indicates we only want to validate input data without creating any new patients and cases.
    """



#        self.run                  = run
#        self.mode                 = mode
#        self.analyses             = analyses
#        if mode == 'germinal':
#            self.termesDF         = self.load_termes(termes)
#        elif mode == 'somatic':
#            self.termesDF         = pd.DataFrame(columns=self.termesDF_column_types.keys()).astype(self.termesDF_column_types)



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
        termesDF = pd.DataFrame(columns=q.termesDF_column_types.keys()).astype(q.termesDF_column_types)
        analysesDF = q.load_analyses(analyses)

        for index, analyse in analysesDF.iterrows():
            print (analyse)
            termesDF = q.reset_analyse_to_termes(analyse, termesDF)

            # Validates the input analyses by iterating members
            family_members = termesDF

            # Iterate trought family members
            for index, individual in family_members.iterrows():
                sample_nanuq = json.loads(nq.get_sample(individual['aliquot']))[0]
                family_members.loc[index] = q.fix_nanuq_termes (sample_nanuq, individual)

            # Create and push analysis
            analysis_payload = q.get_analysis_payload(family_members, args.fix)
            if args.validate: 
                analysis_payload = q.create_analysis(analysis_payload, True)
                print (f"Validation succesfull")
            else: 
                analysis_payload = q.create_analysis(analysis_payload, False)
                # Push sequencings data
                sequencing_payload = q.get_somatic_sequencings_payload(analysis_payload)
                sequencing_response = q.create_sequencings(sequencing_payload)
                # Start pipeline
                pipeline_payload = q.get_pipeline_payload(sequencing_payload)
                pipeline_response = q.create_pipeline(pipeline_payload)

    print (f"Completed")

if __name__ == '__main__':
    main()

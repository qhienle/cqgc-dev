#!/usr/bin/env python3

"""
Utilities for the Nannuq API. Import as a package from another script:

```
from nanuq import Nanuq
nq = Nanuq()
print(nq.get_sample(21057))
```

Or as satndalone script. 

USAGE: 
    nanuq.py --help
    # Download 'SampleSheet.csv', 'SampleNames.txt' and 'SamplePools.csv'
    # for run A00516_0428
    nanuq.py A00516_0428
    # Alternative format for Run name 
    nanuq.py Seq_S2_PRAG_20230811

See also the SampleSheet Class from `samplesheet.py`

## Sample code

Examples from Nanuq:

```bash
wget --post-data "j_username=USER&j_password=PASS" --no-cookies https://cigcp-nanuq.calculquebec.ca/nanuqMPS/sampleSheet/NovaSeq/A00516_0295/ -O 'SampleSheet.csv'
wget --post-data "j_username=${2}&j_password=${3}" --no-cookies "${server}/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
Get credentials from file
echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq
wget --post-file ~/.nanuq --no-cookies "https://nanuq.cqgc.hsj.rtss.qc.ca/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
```

Converted to python, with subprocess. Here, we use 'requests' 

```python
import subprocess
url = https://nanuq.cqgc.hsj.rtss.qc.ca/nanuqMPS/ws/GetClinicalSampleInfoWS?name=21057
subprocess.run(['wget', '--post-file', '~/.nanuq', '--no-cookies', url, '-O', 'nanuq.json'])
```
"""
import os, sys
import argparse
import requests
import datetime
import logging

__version__ = "0.3"

CONFIG_FILE = os.path.expanduser('~') + os.sep + '.nanuq'

    
def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Utilities for Nanuq", epilog="Ex.: python nanuq.py -r A00516_0428")
    parser.add_argument('-r', '--run',      help="Get Nanuq files for Run ID, ex: 'A00516_0428'")
    parser.add_argument('-n', '--names',    help="Get sample names in Run ID, ex: 'A00516_0428'")
    parser.add_argument('-s', '--sample',   help="Get information for Nanuq sample (CQGC) ID, ex: '40890'")
    parser.add_argument('-u', '--username', help="Nanuq username")
    parser.add_argument('-p', '--password', help="Nanuq password")
    parser.add_argument('-n', '--no-check-run-name', action='store_true', dest='skip_check', help="Do not check run name")
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


class Nanuq:
    """
    Return an object to interact with Nanuq
    """
    def __init__(self, username=None, password=None, config_file=CONFIG_FILE):
        """
        Configure settings
        """
        self.config_file = config_file
        self.server      = 'https://nanuq.cqgc.hsj.rtss.qc.ca' # 'http://spxp-app07'
        self.username    = username
        self.password    = password
        self.auth_data   = self.get_auth()


    def configure(self, username, password, config_file=CONFIG_FILE):
        """
        Write username and password to config_file
        """
        self.username    = username
        self.password    = password
        self.config_file = config_file
        with open(config_file, 'w') as cfg:
            cfg.write(f"j_username={username}&j_password={password}&toto=1")
        print(f'Wrote username and password to configuration file {config_file}')


    def get_auth(self):
        """
        Generate authentication data for Nanuq. If username and password info are
        not provided, look for it in configuration file `~/.nanuq`, which contains
        the single line: j_username=USERNAME&j_password=PASSWORD&toto=1
        """
        if self.username is None and self.password is None:
            # Credentials not provided, get them from ~/.nanuq file
            #
            with open(os.path.expanduser("~/.nanuq"), 'r') as auth_fh:
                auth_parts = auth_fh.readline().split('&')
                self.username = auth_parts[0].split('=')[1]
                self.password = auth_parts[1].split('=')[1]
        return {'j_username': self.username, 'j_password': self.password}
    

    def get_api(self, url, outfile=None):
        """
        GET {url} from Nanuq. Returns a requests' response object, and a file,
        if filename for {outfile} is not 'None'.
        """
        try:
            logging.debug(f"Connecting to {url}")
            response = requests.post(url, data=self.auth_data)
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            logging.warning(f"{err}") # print(f"HTTP status code: {response.status_code}")
            if response.status_code == 400:
                logging.warning(f"Bad request. Try using an alternative run format 'A00516_0447' or 'Seq_S2_PRAG_20230811'")
        if outfile is None:
            return(response)
        else:
            with open(outfile, 'w') as samplesheet_fh:
                samplesheet_fh.write(response.text)
            return(response)
    
    
    def check_run_name(self, run):
        """
        Check run name format. Ex.: "A00516_0106" or "Seq_S2_PRAG_20230811".
        If necessary, convert run identifier from long to short (ex: from 
        "200302_A00516_0106_BHNKHFDMXX" to "A00516_0106").
        - run: RunID (ex: "200302_A00516_0106_BHNKHFDMXX" or "A00516_0106")
        - Returns: valid short format for RunID (ex: "A00516_0106")
        """
        sequencers = ['LH00336', 'A00516', 'A00977', 'NB551410', 'LH00207R']
        fc_parts = run.split('_')
        if len(fc_parts) == 4 and fc_parts[0] == 'Seq':
            # Format "Seq_S2_PRAG_20230811"
            #
            logging.warning(f"Run name {run} appears to be experiment name.")
            return(run)
        else:
            # Format "200302_A00516_0106_BHNKHFDMXX" or "A00516_0106".
            # Long form needs to be shortened.
            #
            if len(fc_parts) == 2 and fc_parts[0] in sequencers and len(fc_parts[1]) == 4:
                return(run)
            elif len(fc_parts) == 4 and fc_parts[1] in sequencers and len(fc_parts[2]) == 4:
                # RunID given in long format, ex.: "200302_A00516_0106_BHNKHFDMXX"
                #
                fc_short = f"{fc_parts[1]}_{fc_parts[2]}"
                logging.info(f"Run name {run} in long format. Converted to short form {fc_short}")
                return(fc_short)
            else:
                raise ValueError(f"Incorrect format for RunID {run}. Please use a value like 'A00516_0106'.")
    

    def parse_run_name(self, run):
        """
        Parse run identifier.
        - run    : [str] Illumina Run ID, ex: 20240510_LH00336_0043_A22K5KMLT3
        - Returns: A tuple (date, instrument_id, run_num, fc_id, fc_short),
                   where date is of type [DateTime] and the rest are [str].
                   e.g., ("2024-05-10", "LH00336", "0043", "A22K5KMLT3", 
                   "LH00336_0043").
        """
        fc_short = self.check_run_name(run) # ex: LH00336_0043
        fc_parts = run.split('_')
        if len(fc_parts) == 2:
            logging.warning(f"Run name {run} appears to be in short format")
            return None, fc_parts[0], fc_parts[1], None, fc_short
        elif len(fc_parts) == 4:
            # Better to convert DateTime based on the instrument ID (fc_parts[1])?
            # NovaSeqX (LH00336) has 8 digits for dates (yyyymmdd), 
            # whereas NovaSeq6000 (A00516, A00977) have 6 (yymmdd).
            #
            if fc_parts[0].isdigit():
                if len(fc_parts[0]) == 8:
                    date = datetime.datetime.strptime(fc_parts[0], '%Y%m%d').strftime('%Y-%m-%d')
                elif len(fc_parts[0]) == 6:
                    date = datetime.datetime.strptime(fc_parts[0], '%y%m%d').strftime('%Y-%m-%d')
            else:
                logging.warning(f"Wrong format for flowcell date '{fc_parts[0]}'")
        else:
            logging.info(f"Wrong format for Run name '{run}'")
            date = fc_parts[0]
        return date, fc_parts[1], fc_parts[2], fc_parts[3], fc_short


    def get_samplesheet(self, run, outfile=None, skip_check=False):
        """
        Download SampleSheet for `run` from Nanuq.
        Writes to `outfile`, if specified.
        Returns a requests response object.
        """
        if skip_check is False:
            run = self.check_run_name(run)
        # Index2 in reverse-complement for NovaSeq6000
        # url = f'{self.server}/nanuqMPS/sampleSheetV2/NovaSeq/{run}'
        # New API decides whether Index2 is rc for NovaSeq6000, or forward for NovaSeqX
        url = f'{self.server}/nanuqMPS/dragenSampleSheet/NovaSeq/{run}/'
        response = self.get_api(url, outfile)
        return(response)


    def get_samplenames(self, run, outfile=None, skip_check=False):
        """
        Download SampleNames for `run` from Nanuq.
        Writes to `outfile`, if specified.
        Returns a requests response object, and a SampleNames.txt file.
        """
        if skip_check is False:
            run = self.check_run_name(run)
        url = f'{self.server}/nanuqMPS/sampleConversionTable/run/{run}/technology/NovaSeq/'
        response = self.get_api(url, outfile)
        return(response)


    def get_samplepools(self, run, outfile=None, skip_check=False):
        """
        Download SamplePools for `run` from Nanuq.
        Writes to `outfile`, if specified.
        Returns a requests response object, and a SamplePools.csv file.
        """
        if skip_check is False:
            run = self.check_run_name(run)
        url = f'{self.server}/nanuqMPS/poolingSampleSheet/run/{run}/technology/NovaSeq/'
        response = self.get_api(url, outfile)
        return(response)


    def get_sample(self, id_cqgc):
        """
        Get sample id_cqgc from Nanuq. Returns a JSON string.
        - id_cqgc: [str] e.g. 32994
        - Returns: [str] in JSON format
        """
        url = f'{self.server}/nanuqMPS/ws/GetSampleInfoWS?name={id_cqgc}'
        response = self.get_api(url)
        return(response.text)
    

    def get_clinical_sample(self, id_cqgc):
        """
        Get sample id_cqgc from Nanuq. Returns a JSON string. 
        - id_cqgc: [str] e.g. 22762
        - Returns: [str] in JSON format
        """
        url = f'{self.server}/nanuqMPS/ws/GetClinicalSampleInfoWS?name={id_cqgc}'
        response = self.get_api(url)
        return(response.text)
    

    def list_samples(self, run, file=None):
        """
        Return a list of tuples (cqgc_id, sample_names) for run. 
        A SampleNames.txt `file` can be provided instead of downloading from Nanuq.
        - run (str) : Run name
        - file (str): Nanuq's "SampleNames.txt" or a one-column list of CQGC IDs.
        """
        samples = []
        if file is None:
            samplenames = self.get_samplenames(run)
            if not samplenames.text.startswith("##20"):
                sys.exit(logging.error(f"Unexpected content for SampleNames. Please verify Nanuq's reponse:\n{samplenames.text}"))
            else:
                logging.info("Retrieved samples conversion table from Nanuq")
                lines = samplenames.text.splitlines()
        else:
            logging.info(f"Using list of samples from file {file} instead of Nanuq")
            with open(file, 'r') as fh:
                lines = fh.readlines()

        for line in lines:
            if not line.startswith('#'):
                samples.append(tuple(line.split("\t")))
        return samples


def main():
    args = parse_args()
    os.chdir(os.getcwd())
    nq = Nanuq()
        
    if args.run:
        # Write Nanuq files for given 'run'
        #
        nq.get_samplesheet(args.run, outfile='SampleSheet.csv', skip_check=args.skip_check)
        nq.get_samplenames(args.run, outfile='SampleNames.txt', skip_check=args.skip_check)
        nq.get_samplepools(args.run, outfile='SamplePools.csv', skip_check=args.skip_check)
    elif args.names:
        # Convert names for samples in given run
        #
        fc_parts = nq.parse_run_name(args.names)
        print(fc_parts)
        for t in nq.list_samples(fc_parts[4]):
            print(t)
    elif args.sample:
        # Get detailed information for 'sample'
        #
        print(nq.get_sample(args.sample))
    else:
        print("\nType `nanuq.py --help` for more on usage.\n")


if __name__ == '__main__':
    main()

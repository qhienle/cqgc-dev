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
import os
import argparse
import requests
import time

__version__ = "0.2"

CONFIG_FILE = os.path.expanduser('~') + os.sep + '.nanuq'

    
def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Utilities for Nanuq", epilog="Ex.: python nanuq.py -r A00516_0428")
    parser.add_argument('-r', '--run', required=True, help="Run ID, ex: 'A00516_0428'")
    parser.add_argument('-u', '--username', help="Nanuq username")
    parser.add_argument('-p', '--password', help="Nanuq password")
    parser.add_argument('-n', '--no-check-run-name', action='store_true', dest='skip_check', help="Do not check run name")
    return(parser.parse_args())


def now():
    """
    Returns a timestamp string, ex: print(f'{now} is the right time to salute')
    """
    return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))


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
            # print(f"{now()} Connecting to {url}")
            response = requests.post(url, data=self.auth_data)
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            print(f"{now()} {err}") # print(f"HTTP status code: {response.status_code}")
            if response.status_code == 400:
                print(f"{now()} Bad request. Try using an alternative run format 'A00516_0447' or 'Seq_S2_PRAG_20230811'")
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
        sequencers = ['LH00336', 'A00516', 'A00977', 'NB551410']
        fc_parts = run.split('_')
        if len(fc_parts) == 4 and fc_parts[0] == 'Seq':
            # Format "Seq_S2_PRAG_20230811"
            #
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
                print(f"{now()} RunID {run} in long format. Converted to short form {fc_short}")
                return(fc_short)
            else:
                raise ValueError(f"Incorrect format for RunID {run}. Please use something like 'A00516_0106' or skip_check with `--no-check-run-name`.")
    

    def get_samplesheet(self, run, outfile=None, skip_check=False):
        """
        Download SampleSheet for `run` from Nanuq.
        Writes to `outfile`, if specified.
        Returns a requests response object.
        """
        if skip_check is False:
            run = self.check_run_name(run)
        url = f'{self.server}/nanuqMPS/sampleSheetV2/NovaSeq/{run}'
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
        TODO Or a JSON onj is better?
        """
        url = f'{self.server}/nanuqMPS/ws/GetClinicalSampleInfoWS?name={id_cqgc}'
        response = self.get_api(url)
        return(response.text)
    

def main():
    args = parse_args()
    os.chdir(os.getcwd())
    nanuq = Nanuq()
    nanuq.get_samplesheet(args.run, outfile='SampleSheet.csv', skip_check=args.skip_check)
    nanuq.get_samplenames(args.run, outfile='SampleNames.txt', skip_check=args.skip_check)
    nanuq.get_samplepools(args.run, outfile='SamplePools.csv', skip_check=args.skip_check)
    print("\nDone.\n")


if __name__ == '__main__':
    main()

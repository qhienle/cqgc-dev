#!/usr/bin/env python3
"""
Get files from Nanuq

Download `SampleSheet.csv`, `SampleNames.txt` and `SamplePools.csv` for a given
_Run_ from _Nanuq_.

USAGE: get_nanuq_files.py
       get_nanuq_files.py [-u USERNAME -p PASSWORD] -r RUN
       get_nanuq_files.py --help

<RUN> has to be in FC_SHORT format, _e.g._ 'A00516_0339'. If not provided, the
shell's global environment variable ${FC_SHORT} is used (defined at the start 
of the procedure for bcl-convert).

Nanuq username and password can be saved to a file named '~/.nanuq':
`echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`
Replace USERNAME and PASSWORD with actual values.

Files are saved to the current working directory, or to `${WORKDIR}/${FC}`, 
if these environment variables are set globally _e.g._:
`export FC="220224_A00516_0331_AHKLTJDRXY"`.
"""

import os, sys, subprocess
import argparse


# TODO: Move functionalities to Nanuq Class, instead of subprocess
# 
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from lib.nanuq import Nanuq


__version__ = "0.1"

def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Download files for a Run from Nanuq")
    parser.add_argument('-r', '--run',      help="FC_SHORT Run ID, ex: 'A00516_339'")
    parser.add_argument('-u', '--username', help="Nanuq username")
    parser.add_argument('-p', '--password', help="Nanuq password")
    return(parser.parse_args())

def download_files(run, credentials, out_sheet, out_names, out_pools):
    """
    Examples from Nanuq
    wget --post-data "j_username=USER&j_password=PASS" --no-cookies https://cigcp-nanuq.calculquebec.ca/nanuqMPS/sampleSheet/NovaSeq/A00516_0295/ -O 'SampleSheet.csv'
    wget --post-data "j_username=${2}&j_password=${3}" --no-cookies "${server}/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
    Get credentials from file
    echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq
    wget --post-file ~/.nanuq --no-cookies "https://nanuq.cqgc.hsj.rtss.qc.ca/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
    """
    #server = 'http://spxp-app07'
    server = 'https://nanuq.cqgc.hsj.rtss.qc.ca'

    url_sheet = f'{server}/nanuqMPS/sampleSheetV2/NovaSeq/{run}/'
    url_names = f'{server}/nanuqMPS/sampleConversionTable/run/{run}/technology/NovaSeq/'
    url_pools = f'{server}/nanuqMPS/poolingSampleSheet/run/{run}/technology/NovaSeq/'

    subprocess.run(['wget', '--post-data', credentials, '--no-cookies', url_sheet, '-O', out_sheet])
    subprocess.run(['wget', '--post-data', credentials, '--no-cookies', url_names, '-O', out_names])
    subprocess.run(['wget', '--post-data', credentials, '--no-cookies', url_pools, '-O', out_pools])

def main(run=None, username=None, password=None):
    """
    Connect to Nanuq web API and download three files.
    If Nanuq credentials provided as params, download files.
    Otherwise, look for them in config file `~/.nanuq`.
    """
    auth = '~/.nanuq'

    # If RUN identifier not defined by option -r/--run, try to get this info
    # from env var ${FC_SHORT}, if set as part of the bcl-convert procedure.
    #
    if run is None:
        print(f"\nWARNING: Option -r/--run not provided. Get FC_SHORT from environment settings.")
        try:
            run = os.environ['FC_SHORT']
        except KeyError:
            raise SystemExit("ERROR: Option -r/--run not provided and ${FC_SHORT} not set.\n")
        except Exception as err:
            raise SystemExit(f"Caught an unexpected Exception:\n{err}\n")
    else:
        print(f"Got FC_SHORT '{run}' from command line argument.")

    # Connect and download, using Nanuq authentication with credentials from 
    # either the command-line arguments or from the config file ./nanuq
    #
    if username is not None or password is not None:
        credentials = f"j_username={username}&j_password={password}"
    else:
        # Username and password not provided. Look for them in ~/.nanuq
        #
        try:
            print(f"\nLooking for credentials in file '{auth}'")
            with open(os.path.expanduser(auth), 'r') as fh:
                line = fh.readline().rstrip()
        except FileNotFoundError as fnf:
            raise SystemExit(f"{fnf}.\nPlease provide -u/--username and -p/--password for Nanuq.\n")
        except Exception as err:
            raise SystemExit(f"Caught an unexpected Exception:\n{err}\n")
        else:
            # Extract credentials from file '~/.nanq', which should contain a
            # single line: 'j_username=USERNAME&j_password=PASSWORD&toto=1'
            # Remove last part, `toto=1`, only needed for `wget --post-file`
            #
            parts = line.split('&')
            credentials = parts[0] + '&' + parts[1]

    # Set output directory. If the shell global environment variable for Run
    # is set, define output folder from those values: ${WORKDIR}/${FLOWCELL}/
    # Otherwise, save files to current working directory.
    #
    if 'WORKDIR' in os.environ and 'FC' in os.environ:
        outdir = os.environ['WORKDIR'] + os.sep + os.environ['FC']
        os.makedirs(outdir, exist_ok=True)
    else:
        print(f"\nWARNING: Environment not set. Downloads in: {os.getcwd()}.\n")
        outdir = os.getcwd() + os.sep

    # Download the three files, with the parameters gathered above.
    #
    out_sheet = outdir + os.sep + 'SampleSheet.csv'
    out_names = outdir + os.sep + 'SampleNames.txt'
    out_pools = outdir + os.sep + 'SamplePools.csv'

    download_files(run, credentials, out_sheet, out_names, out_pools)

    # If files downloaded files are empty (size == 0), no run could be found 
    # with ${FC_SHORT} identifier. Try with ${XP}, if set, else warn and quit. 
    #
    if os.stat(out_sheet).st_size == 0 and os.stat(out_names).st_size == 0 and os.stat(out_pools).st_size == 0:
        print(f"ERROR: Empty files downloaded with {run}. Trying with ${{XP}}.\n")
        if 'XP' in os.environ:
            run = os.environ['XP']
            download_files(run, credentials, out_sheet, out_names, out_pools)
        else:
            print(f"ERROR: Empty files downloaded with {run} and ${{XP}} not set.")
            print(f"Please check identifier for RUN.")

    print(f"First lines of downloaded files:")
    #subprocess.run(['head', out_sheet, out_names, out_pools])
    for file in [out_sheet, out_names, out_pools]:
        print(f"\n===> {file}:\n")
        with open(file, 'r', encoding='ISO-8859-1') as fh:
            for line in fh.readlines()[:5]:
                print(line.rstrip())

if __name__ == '__main__':
    args = parse_args()
    main(args.run, args.username, args.password)
    print("\nDone.\n")

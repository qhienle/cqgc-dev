#!/usr/bin/env python3
"""
Get files from Nanuq

Download from Nanuq `SampleSheet.csv`, `SampleNames.txt` and `SamplePools.csv` 
for a given Run (a.k.a FLOWCELL).

USAGE: get_nanuq_files.py [-u <USERNAME> -p <PASSWORD>] -r <RUN>
       get_nanuq_files.py -r 20250314_LH00336_0185_B22WGWFLT3
       get_nanuq_files.py --help

Nanuq username and password can be saved to a file named '~/.nanuq':
`echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`
Replace USERNAME and PASSWORD with actual values.

Files are saved to the current working directory, or to `${WORKDIR}/${FC}`, 
if these environment variables are set globally _e.g._:
`export FC="220224_A00516_0331_AHKLTJDRXY"`.

## OPTIONS

--orient-index2|-o: [DEPRECATED] Orient index2 in the SampleSheet for NovaSeqX
    (sense). If FALSE (not present), index2 in the SampleSheet will be in 
    reverse-complement (for NovaSeq6000). 
    Orientation of index2 is now determined by the instrument identifier in the
    RUN name (e.g. A00516=NovaSeq6000, LH00336=NovaSeqX)
"""

import os, sys, subprocess
import argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from lib.nanuq import Nanuq


__version__ = "0.1"

def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Download from Nanuq files for a Run/Flowcell")
    parser.add_argument('-r', '--run',      help="Run ID, ex: '250115_A00516_0640_BHLJG3DSXC'")
    parser.add_argument('-u', '--username', help="Nanuq username")
    parser.add_argument('-p', '--password', help="Nanuq password")
    parser.add_argument('-n', '--no-check-run-name', dest='no_check', action="store_true", 
                        help="Do not verify validity of -r/--run name")
    parser.add_argument('-o', '--orient-index2', action="store_true", 
                        help="(DEPRECATED) Orient index2 in the SampleSheet for NovaSeqX")
    return(parser.parse_args())

def download_files(run, credentials, out_sheet, out_names, out_pools, no_check=False):
    """
    Download SampleSheet.csv, SampleNames.txt and SamplePools.csv from Nanuq
    - run [str]: short Run name, ex: LH00336_0043
    - credentials [str]: username&password. See below
    - out_sheet: SampleSheet file. Default=SampleSheet.csv
    - out_names: SampleNames file. Default=SampleNames.txt
    - out_sheet: SamplePools file. Default=SamplePools.csv
    - no_check : Do not verify validity of run name
        
    Examples from Nanuq
    wget --post-data "j_username=USER&j_password=PASS" --no-cookies https://cigcp-nanuq.calculquebec.ca/nanuqMPS/sampleSheet/NovaSeq/A00516_0295/ -O 'SampleSheet.csv'
    wget --post-data "j_username=${2}&j_password=${3}" --no-cookies "${server}/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
    Get credentials from file
    echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq
    wget --post-file ~/.nanuq --no-cookies "https://nanuq.cqgc.hsj.rtss.qc.ca/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
    """
    server = 'https://nanuq.cqgc.hsj.rtss.qc.ca' # 'http://spxp-app07'
    nq = Nanuq()
    if no_check:
        fc_short = run
    else:
        fc_short = nq.check_run_name(run)
    instrument = fc_short.split('_')[0]

    # Different Nanuq API endpoints determine whether Index2 is in reverse-
    # complement for NovaSeq6000 (A00 series of instrument IDs), or forward 
    # for NovaSeqX (LH00 series)
    #
    if instrument.startswith('A00'):
        url_sheet = f'{server}/nanuqMPS/sampleSheetV2/NovaSeq/{fc_short}/'
    elif instrument.startswith('LH00'):
        url_sheet = f'{server}/nanuqMPS/dragenSampleSheet/NovaSeq/{fc_short}/'
    else:
        sys.exit(f"ERROR: Could not determine index2 orientation for {instrument}")
    url_names = f'{server}/nanuqMPS/sampleConversionTable/run/{fc_short}/technology/NovaSeq/'
    url_pools = f'{server}/nanuqMPS/poolingSampleSheet/run/{fc_short}/technology/NovaSeq/'

    subprocess.run(['wget', '--post-data', credentials, '--no-cookies', url_sheet, '-O', out_sheet])
    subprocess.run(['wget', '--post-data', credentials, '--no-cookies', url_names, '-O', out_names])
    subprocess.run(['wget', '--post-data', credentials, '--no-cookies', url_pools, '-O', out_pools])

def main():
    """
    Connect to Nanuq web API and download three files.
    If Nanuq credentials provided as params, download files.
    Otherwise, look for them in config file `~/.nanuq`.
    """
    args = parse_args()
    auth = f"{os.path.expanduser('~')}{os.sep}.nanuq"

    # If RUN identifier not defined by option -r/--run, try to get this info
    # from env var ${FC_SHORT}, if set as part of the bcl-convert procedure.
    #
    if args.run is None:
        print(f"\nWARNING: Option -r/--run not provided. Get Run ID (FC) from environment settings.")
        try:
            fc = os.environ['FC']
        except KeyError:
            raise SystemExit("ERROR: Option -r/--run not provided and ${FC} not set.\n")
        except Exception as err:
            raise SystemExit(f"Caught an unexpected Exception:\n{err}\n")
    else:
        fc = args.run
        print(f"Got Run ID (FC) '{fc}' from command line argument.")

    # WARN that CLI option `--orient-index2` is DEPRECATED
    #
    if args.orient_index2:
        print(f"WARNING: Option --orient-index2 is DEPRECATED. Orientation of index2 determined from RUN name")

    # Connect and download, using Nanuq authentication with credentials from 
    # either the command-line arguments or from the config file ./nanuq
    #
    if args.username is not None or args.password is not None:
        credentials = f"j_username={args.username}&j_password={args.password}"
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

    download_files(fc, credentials, out_sheet, out_names, out_pools, no_check=args.no_check)

    # If files downloaded files are empty (size == 0), no run could be found 
    # with ${FC_SHORT} identifier. Try with ${XP}, if set, else warn and quit. 
    #
    if os.stat(out_sheet).st_size == 0 and os.stat(out_names).st_size == 0 and os.stat(out_pools).st_size == 0:
        print(f"ERROR: Empty files downloaded with {fc}. Trying with ${{XP}}.\n")
        if 'XP' in os.environ:
            download_files(os.environ['XP'], credentials, out_sheet, out_names, out_pools)
        else:
            print(f"ERROR: Empty files downloaded with {fc} and ${{XP}} not set.")
            print(f"Please check identifier for RUN.")

    print(f"First lines of downloaded files:")
    for file in [out_sheet, out_names, out_pools]:
        print(f"\n===> {file}:\n")
        with open(file, 'r', encoding='ISO-8859-1') as fh:
            for line in fh.readlines()[:5]:
                print(line.rstrip())

if __name__ == '__main__':
    main()
    print("\nDone.\n")

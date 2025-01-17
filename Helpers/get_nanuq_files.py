#!/usr/bin/env python3
"""
Get files from Nanuq

Download `SampleSheet.csv`, `SampleNames.txt` and `SamplePools.csv` for a given
_Run_ from _Nanuq_.

USAGE: get_nanuq_files.py -r LH00336_0043
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

## OPTIONS

--orient-index2|-o: [DEPRECATED] Orient index2 in the SampleSheet for NovaSeqX
    (sense). If FALSE (not present), index2 in the SampleSheet will be in 
    reverse-complement (for the NovaSeq6000). 
"""

import os, sys, subprocess
import argparse

#sys.path.append('/staging2/soft/CQGC-utils')
#src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
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
    parser.add_argument('-o', '--orient-index2', action="store_true", 
                        help="Orient index2 in the SampleSheet for NovaSeqX (deprecated)")
    return(parser.parse_args())

def download_files(run, credentials, out_sheet, out_names, out_pools, orient_index2=False):
    """
    Download SampleSheet.csv, SampleNames.txt and SamplePools.csv from Nanuq
    - run [str]: short Run name, ex: LH00336_0043
    - credentials [str]: username&password. See below
    - out_sheet: SampleSheet file. Default=SampleSheet.csv
    - out_names: SampleNames file. Default=SampleNames.txt
    - out_sheet: SamplePools file. Default=SamplePools.csv
    - orient_index2: API call to select the version of the SampleSheet, with 
        index2 in sense (NovaSeqX) or reverse-complement (NovaSeq600).
        Default=False (reverse-complement, for the NovaSeq6000)
        
    Examples from Nanuq
    wget --post-data "j_username=USER&j_password=PASS" --no-cookies https://cigcp-nanuq.calculquebec.ca/nanuqMPS/sampleSheet/NovaSeq/A00516_0295/ -O 'SampleSheet.csv'
    wget --post-data "j_username=${2}&j_password=${3}" --no-cookies "${server}/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
    Get credentials from file
    echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq
    wget --post-file ~/.nanuq --no-cookies "https://nanuq.cqgc.hsj.rtss.qc.ca/nanuqMPS/sampleSheetV2/NovaSeq/A00516_0339/" -O "SampleSheet.csv"
    """
    #server = 'http://spxp-app07'
    server = 'https://nanuq.cqgc.hsj.rtss.qc.ca'
    nq = Nanuq()

    # Different Nanuq API endpoints determine whether Index2 is in reverse-
    # complement for NovaSeq6000 (A00 series of instrument IDs), or forward 
    # for NovaSeqX (LH00 series)
    #
    instrument = nq.parse_run_name(run)[1]
    if instrument.startswith('A00'):
        url_sheet = f'{server}/nanuqMPS/sampleSheetV2/NovaSeq/{run}/'
    elif instrument.startswith('LH00'):
        url_sheet = f'{server}/nanuqMPS/dragenSampleSheet/NovaSeq/{run}/'
    else:
        sys.exit(f"ERROR: Could not determine index2 orientation for {instrument}")
    url_names = f'{server}/nanuqMPS/sampleConversionTable/run/{run}/technology/NovaSeq/'
    url_pools = f'{server}/nanuqMPS/poolingSampleSheet/run/{run}/technology/NovaSeq/'

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
    #auth = '~/.nanuq'
    auth = f"{os.path.expanduser('~')}{os.sep}.nanuq"

    # If RUN identifier not defined by option -r/--run, try to get this info
    # from env var ${FC_SHORT}, if set as part of the bcl-convert procedure.
    #
    if args.run is None:
        print(f"\nWARNING: Option -r/--run not provided. Get FC_SHORT from environment settings.")
        try:
            args.run = os.environ['FC_SHORT']
        except KeyError:
            raise SystemExit("ERROR: Option -r/--run not provided and ${FC_SHORT} not set.\n")
        except Exception as err:
            raise SystemExit(f"Caught an unexpected Exception:\n{err}\n")
    else:
        print(f"Got FC_SHORT '{args.run}' from command line argument.")

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

    print(f"\n\nCREDENTIALS: {credentials}; OUTDIR: {outdir}; OUT_SHEET: {out_sheet}\n\n")
    download_files(args.run, credentials, out_sheet, out_names, out_pools, args.orient_index2)

    # If files downloaded files are empty (size == 0), no run could be found 
    # with ${FC_SHORT} identifier. Try with ${XP}, if set, else warn and quit. 
    #
    if os.stat(out_sheet).st_size == 0 and os.stat(out_names).st_size == 0 and os.stat(out_pools).st_size == 0:
        print(f"ERROR: Empty files downloaded with {args.run}. Trying with ${{XP}}.\n")
        if 'XP' in os.environ:
            download_files(os.environ['XP'], credentials, out_sheet, out_names, out_pools, args.orient_index2)
        else:
            print(f"ERROR: Empty files downloaded with {args.run} and ${{XP}} not set.")
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

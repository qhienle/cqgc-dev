#!/usr/bin/env python3
"""
Make an "analysis.txt" file from sample information in Nanuq for a given run.

USAGE: make_analysis_file.py RUN_ID [-a], [-u -p], [-i]
       make_analysis_file.py --help
       ex: make_analysis_file.py 230601_A00516_0422_BH5T55DRX3
           make_analysis_file.py 230601_A00516_0422_BH5T55DRX3 -u USERNAME -p PASSWORD
           make_analysis_file.py 230601_A00516_0422_BH5T55DRX3 -i SampleNames.txt -a WES_germinal

Nanuq authentication is expected to be found in a file named '~/.nanuq'.
`echo "j_username=USERNAME&j_password=PASSWORD&toto=1" > ~/.nanuq`

Sample information downloaded from Nanuq is saved to a file 'SampleNames.txt'.
SampleNames can also be provided from a file instead of being fetched using the
`-f|--file` optional argument (default='./SampleNames.txt').

See examples of SampleNames.txt (input) and analysis.txt (output) files in the
annex at the end of this script.
"""

import os
import argparse
import time
import subprocess
import requests

__version__ = "0.1"


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Make an 'analysis.txt' file from Nanuq information")
    parser.add_argument('run', 
                        help="Run ID (ex: 230601_A00516_0422_BH5T55DRX3) [REQUIRED]")
    parser.add_argument('-a', '--analysis-type', 
                        dest='analysis',
                        choices=['RNASeq_somatic', 'WES_germinal', 'WES_somatic-tumor_only', 'WES_somatic-tumor_normal'],
                        default='RNASeq_somatic',
                        help="Name of analysis type (pipeline to use). Default=RNASeq_somatic")
    parser.add_argument('-i', '--input-file', 
                        dest='infile',
                        help="Get information from file instead of download from Nanuq. Default=SampleNames.txt")
    parser.add_argument('-u', '--username', 
                        help="Username for Nanuq")
    parser.add_argument('-p', '--password', 
                        help="Password for Nanuq")
    parser.add_argument('-t', '--test', action="store_true", help="Run _test function only")

    return(parser.parse_args())


def now():
    """
    Returns a timestamp string, ex: print(f'{now} is the right time to salute')
    """
    return(time.strftime('[%Y-%m-%d@%H:%M:%S] '))


def auth_nanuq(username, password):
    """
    Generate authentication data for Nanuq. If username and password info are
    not provided, look for it in configuration file `~/.nanuq`, which contains
    the single line: j_username=USERNAME&j_password=PASSWORD&toto=1
    """
    if username is None and password is None:
        # Credentials not provided, get them from ~/.nanuq file
        #
        with open(os.path.expanduser("~/.nanuq"), 'r') as auth_fh:
            auth_parts = auth_fh.readline().split('&')
            username = auth_parts[0].split('=')[1]
            password = auth_parts[1].split('=')[1]
    return {'j_username': username, 'j_password': password}


def convert_nanuq2analysis(infile, run, analysis):
    """
    Convert Nanuq's 'SampleNames.txt' file to 'analysis.txt'
    - infile  : input file name, e.g.: SampleNames.txt
    - run     : Run ID, ex: 230601_A00516_0422_BH5T55DRX3
    - analysis: type of analysis
    - retuns  : Name of output file (analysis.txt) or None.
    """
    # TODO: try/except opening files
    outfile = os.getcwd() + os.sep + 'analyses.txt'
    with open(infile, 'r') as samplenames_fh, open(outfile, 'w') as analyses_fh:
        for line in samplenames_fh.readlines():
            line = line.rstrip('\r\n')
            # Metadata lines are preceded by '##'
            #
            if line.startswith('##'):
                # Grab the Date
                #
                if line.startswith('##20'):
                    date = line.replace('##', '')
                elif line.startswith('##Flow Cell: '):
                    fc = line.replace('##Flow Cell: ', '')
            else:
                lab_id, sample_id = line.split()
                analyses_str = '\t'.join([date, run, lab_id, 'NA', 'NA', sample_id, sample_id, f"{sample_id}.{analysis}", analysis])
                print(now() + analyses_str)
                analyses_fh.write(analyses_str + '\n')
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        return outfile
    else:
        return None


def main(run, infile, analysis, username=None, password=None):
    """
    Grab the Run ID
    Check analysis type
    Parse the SampleNames file
    Construct the output
    - infile  : input file name, e.g.: SampleNames.txt
    - run     : Run ID, ex: 230601_A00516_0422_BH5T55DRX3
    - analysis: type of analysis
    - retuns  : None. Writes an analysis.txt file.
    """

    workdir = os.getcwd()

    # Short name for run ID is needed to fecth from Nanuq
    # Ex: run=230601_A00516_0422_BH5T55DRX3, fc_short=A00516_0422
    #
    fc_parts = run.split('_')
    fc_short = fc_parts[1] + '_' + fc_parts[2]

    # Analysis type must be one of: 'RNASeq_somatic'
    #
    if analysis is None:
        analysis = 'RNASeq_somatic'
    else:
        pass

    if infile is not None:
        print(f'{now()} Reading samples information from {infile}')
        outfile = convert_nanuq2analysis(infile, run, analysis)
    else:
        print(f'{now()} Downloading SampleNames from Nanuq for {fc_short}')

        # Get a Nanuq authentication string, assemble the request and 
        # download SampleNames.txt from Nanuq
        #
        auth     = auth_nanuq(username, password)
        server   = "https://nanuq.cqgc.hsj.rtss.qc.ca"
        endpoint = f"/nanuqMPS/sampleConversionTable/run/{fc_short}/technology/NovaSeq/"
        url      = server + endpoint
        tmpfile  = workdir + os.sep + 'SampleNames.tmp'

        #subprocess.run(['wget', '--post-data', auth, '--no-cookies', url, '-O', 'SampleNames.tmp'])
        try:
            print(f"{now()} Connecting to {url}")
            response = requests.post(url, data=auth)
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            print(f"{now()} Oops, something went wrong: HTTP status code: {response.status_code}")
            raise SystemExit(err)
        else:
            with open(tmpfile, 'w') as tmp_fh:
                tmp_fh.write(response.text)
            outfile = convert_nanuq2analysis(tmpfile, run, analysis)
            os.remove(tmpfile)

    if outfile is not None:
        print(f"{now()} Wrote {outfile}")
    else:
        print(f"`{now()} Oops, something went wrong: could not create analysis output file.")
    
    return None


def _test(args, opt="."):
    print(f"Required command-line argument is: {args}")
    opt = os.getcwd()
    print(os.stat(opt))
    return None


if __name__ == '__main__':
    args = parse_args()
    if args.test:
        _test(args)
    else:
        main(args.run, args.infile, args.analysis, args.username, args.password)
    print("\nDone.\n")


"""
Example of a SampleNames.txt file (input)

##2023-06-01
##Centre for Pediatric Clinical Genomics
##Flow Cell: H5T55DRX3
##Principal Investigator: Dr Jean-François Soucy
##Nanuq References: A00516_0422
##Content: Internal_Sample_ID -> Client_Sample_Name Conversion grid
##-------------------------------------------
##Internal_Sample_ID	Client_Sample_Name
##-------------------------------------------
21633	OM230223
21634	OM230496
21635	OM230499
21636	OM230551
21637	OM230585
21638	OM160553
##-------------------------------------------
##Description of the conversion grid:
##A "Client_Sample_Name" represents the name of a sequencing sample initially assigned by the client.
##The "Internal_Sample_ID" is the internal identifier assigned by the Center for this sample.
##Use of an "Internal_Sample_ID" ensures the traceability of this sample throughout the sequencing
##process, but also the anonymization of the information generated and transferred in the
##frame of this project.

Example of an analysis.txt file (output)

6/2/2023        230601_A00516_0422_BH5T55DRX3   21633   NA      NA      OM230223        OM230223        OM230396.RNASeq_somatic RNASeq_somatic
6/2/2023        230601_A00516_0422_BH5T55DRX4   21634   NA      NA      OM230496        OM230496        OM230442.RNASeq_somatic RNASeq_somatic
6/2/2023        230601_A00516_0422_BH5T55DRX5   21635   NA      NA      OM230499        OM230499        OM230447.RNASeq_somatic RNASeq_somatic
6/2/2023        230601_A00516_0422_BH5T55DRX6   21636   NA      NA      OM230551        OM230551        OM230448.RNASeq_somatic RNASeq_somatic
6/2/2023        230601_A00516_0422_BH5T55DRX7   21637   NA      NA      OM230585        OM230585        OM230457.RNASeq_somatic RNASeq_somatic
6/2/2023        230601_A00516_0422_BH5T55DRX8   21638   NA      NA      OM160553        OM160553        OM230533.RNASeq_somatic RNASeq_somatic
"""
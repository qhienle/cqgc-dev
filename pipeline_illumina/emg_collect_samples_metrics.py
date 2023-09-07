#!/usr/bin/env python3
"""
Get metrics for PRAGMatIQ samples.

USAGE: emg_get_samples_metrics.py --help

Parse Emedgene's logs files and folders to get analysis metrics for PRAGMatIQ
samples. Log files are located here: 

narval.calculquebec.ca:~/projects/ctb-rallard/COMMUN/PRAGMatIQ-EMG

Logs must first be downloaded from `aws` with the script `archive_PRAGMatIQ.sh`
or using the following command:

`aws s3 --profile emedgene cp s3://cac1-prodca-emg-auto-results/CHU_Sainte_Justine/${sample}/ ./${sample} --recursive`
"""

import os
import argparse
import time

__version__ = 0.1


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Get metrics for PRAGMatIQ samples.")
    parser.add_argument('arg', help="Mandatory argument [REQUIRED]")
    return(parser.parse_args())


def now():
    """
    Returns a timestamp string, ex: print(f"{now} Time to say Hello World!")
    """
    return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))


def main(args):
    """
    From a list of samples, retrieve the following metrics:
    - Institution
    - CQGC_ID
    - Sample_ID
    - #Reads
    - Average Coverage
    - CoverageUniformity
    - %Bases With Coverage > 20X
    - #SNV
    - #Indels
    - #CNVs
    """
    print(f"Required command-line argument is: {args}")
    print("\nDone.\n")


def _test(arg, opt="."):
    print(f"Required command-line argument is: {arg}")
    return(os.stat(opt))


if __name__ == '__main__':
    args = parse_args()
    main(args)
    _test(args)

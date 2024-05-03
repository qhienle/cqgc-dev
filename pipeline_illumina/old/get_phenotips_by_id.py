#!/usr/bin/env python3
"""
Get HPO annotations for Phenotips records at https://chusj.phenotips.com/

USAGE: get_phenotips_by_id.py P0000789
       get_phenotips_by_id.py --help

Get HPO annotations for a given Phenotips record (PID, _e.g.: P0000789).
Prints list of phenotypes identifiers and text labels to STDOUT.

Requires that credentials are stored in a configuration file.
"""

__version__ = "0.2"

import os
import sys
import argparse

# Set source path to CQGC-utils so that we can use relative imports
#
src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(src_path)
from lib.gapp import Phenotips


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Get record from Phenotips by ID, ex: 'P0000789'")
    parser.add_argument('phenotips_id', help="Phenotips identifier, ex: 'P0000789' [REQUIRED]")
    parser.add_argument('-c', '--config-file', dest='config_file', nargs="?", 
                        default=os.path.expanduser("~/.illumina/gapp_conf.json"), 
                        help="Configuration file. Default='~/.illumina/gapp_conf.json'")
    parser.add_argument('-t', '--test', action='store_true', help="Run test only")
    return(parser.parse_args())


def main(args):
    """
    Get patient record by Phenotips ID and prints annotations to STDOUT.
    """
    pho    = Phenotips(config_file=args.config_file)
    try:
        pids = pho.get_hpo(args.phenotips_id)
    except Exception as e:
        raise SystemExit(f"Could not retrieve {args.phenotips_id}: {e}")
    ids    = []
    labels = []
    for pid in pids:
        labels.append(pid['label'])
        ids.append(pid['id'])
        print(pid)
    print("\n", ', '.join(ids))
    print("\n", ', '.join(labels))


def tests(args):
    print("Running tests")
    return(args)


if __name__ == '__main__':
    args = parse_args()
    if args.test:
        tests(args)
    else:
        main(args)
    print(f"\nDone.\n{__file__} version {__version__}\n")

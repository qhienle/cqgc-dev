#!/usr/bin/env python3
"""
Get HPO terms from Phenotips for PID

Retrieve the HPO identifiers and terms (labels) for a Phenotips PID.

USAGE: get_phenotips_hpos.py PID
       get_phenotips_hpos.py --help

<PID> is a Phenotips identifier, _e.g._: P0000808

Phenotips credentials can be saved to a file named '~/.illumina/gapp_conf.json':
{
    "instance"         : "cac1.trusight.illumina.com",
    "X-ILMN-Domain"    : "chusj",
    "X-ILMN-Workgroup" : "blah-blah-blah",
    "X-Auth-Token"     : "APIKey blahblahblah",
    "testDefinitionId" : "blah-blah-blah",
    "bs_apiServer"     : "https://api.cac1.sh.basespace.illumina.com",
    "bs_accessToken"   :  "blah-blah-blah",
    "X-Gene42-Server"  : "https://chusj.phenotips.com",
    "X-Gene42-Auth"    : "Basic tokenblahblahblah",
    "X-Gene42-Secret"  : "secretblahblahblah"
}
"""

import os, sys
import argparse


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from lib.gapp import Phenotips


__version__ = "0.1"

def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Get HPO terms from Phenotips for PID")
    parser.add_argument('pid', help="Phenotips ID, ex.: P0000808")
    return(parser.parse_args())


def main(pid):
    """
    Main stuff
    """


if __name__ == '__main__':
    args = parse_args()
    main(args.pid)
    print("\nDone.\n")

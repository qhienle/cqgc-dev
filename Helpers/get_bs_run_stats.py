#!/usr/bin/env python3
"""
Get statistics for Runs from BaseSpace

Template for Python developments.

USAGE: get_bs_run_stats.py --help
"""

import os
import argparse
import subprocess
import json

__version__ = "0.1"


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Get statistics for Runs from BaseSpace")
    parser.add_argument('run', help="Run name or 'all', ex: 20240911_LH00336_0096_B22NJCNLT3 [str, REQUIRED]")
    return(parser.parse_args())


def bs_list_runs():
    """
    Get statistics from BaseSpace using `bs`
    - arguments:
    - returns:
    """
    run_json = json.loads(
        subprocess.run(['bs', '-c', 'cac1', 'run', 'list', '--format', 'json', '--newer-than', '1d'], 
                       capture_output=True, 
                       text=True
                       ).stdout)
    return run_json


def main(args):
    """
    Define Function1, _e.g._ for testing stuff here
    - arguments:
    - returns:
    """
    bs_list_runs()
    return 1


if __name__ == '__main__':
    main(parse_args())
    print("\nDone.\n")
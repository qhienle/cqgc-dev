#!/usr/bin/env python3
"""
Get statistics for Runs from BaseSpace

Template for Python developments.

USAGE: get_bs_run_stats.py --help
"""

import os
import argparse
import logging
import subprocess
import json

__version__ = "0.1"


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Get statistics for Runs from BaseSpace")
    parser.add_argument('run', help="Run name or 'all', ex: 20240911_LH00336_0096_B22NJCNLT3 [str, REQUIRED]")
    parser.add_argument('--newer-than', '-n', dest='newer', 
                        help="[str] Filter for items that are newer than the \
                        given duration, ex: 5m (5 minutes), 12h. Permitted \
                        suffixes are s, m, h, d, w, y. Assumes that run='all'")
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


def bs_list_runs(run, newer=None):
    """
    Get statistics from BaseSpace using `bs`
    - run:     [str] Run name, or all
    - newer:   [str] Filter items that are newer than [str]. See parse_args().
    - returns: [list] List of runs, where each run is a dict
    """
    if run == 'all':
        if newer is None:
            bs = subprocess.run(['bs', '-c', 'cac1', 'run', 'list', '--format', 'json'], 
                                capture_output=True, text=True).stdout
        else:
            bs = subprocess.run(['bs', '-c', 'cac1', 'run', 'list', '--format', 'json', '--newer-than', newer], 
                                capture_output=True, text=True).stdout
    else:
        bs = subprocess.run(['bs', '-c', 'cac1', 'run', 'list', '--format', 'json', '--filter-term', run],
                            capture_output=True, text=True).stdout
    runs = json.loads(bs)
    print(f"bs command returned {runs}")    
    
    # `bs` returns a run as a dict. If more than one is found, we get a list of
    # dicts. For consistency, we always return a list of zero or more dicts.
    #
    if runs is list:
        return runs
    elif runs is dict:
        return [runs]
    else:
        return []


def main(args):
    """
    Define Function1, _e.g._ for testing stuff here
    - arguments:
    - returns:
    """
    logging.info(f"Listing runs from BaseSpace {args}")
    print(f"Value for newer is {args.newer}.")
    runs = bs_list_runs(run=args.run, newer=args.newer)
    for run in runs:
        pass
    

if __name__ == '__main__':
    main(parse_args())
    logging.info("Done.\n")
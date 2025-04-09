#!/usr/bin/env python3

"""
assemble_report.py

Assemble DRAGEN CSV Reports from `bcl-convert` into sheets of an Excel file.

USAGE: assemble_report.py <Reports_dir>
       assemble_report.py ./210602_A00516_0217_FOO7TESTXX/Reports
       assemble_report.py --help
"""

import sys
import os
import argparse
import logging

__version__ = "0.2"

class Class:
    """
    Define Class here
    """
    def __init__(self):
        """
        Initialize Class
        """
        pass


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Assemble DRAGEN CSV Reports into sheets of an Excel file.")
    parser.add_argument("reports_dir", help="Reports folder with CSV files")
    parser.add_argument("-p", "--pools", help="File containing the pooling fractions. Default=`SamplePools.csv`.")
    parser.add_argument('--logging-level', '-l', dest='level', default='info',
                        help="Logging level (str), can be 'debug', 'info', 'warning'. Default='info'")
    return parser.parse_args()


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


def main(args):
    """
    Main function
    """
    args = parse_args()
    configure_logging(args.level)


if __name__ == '__main__':
    sys.exit(main())

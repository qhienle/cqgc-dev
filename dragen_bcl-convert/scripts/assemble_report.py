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


def main():
    """
    Write a `Reports.xlsx` file and generate a PNG bar plot of observed versus
    expected number of reads. Excel file is a collection of sheets, gathered
    from reports generated by DRAGEN (`<DRAGEN_OUTPUTDIR>/Reports/*.csv`).
    Requires the default "Reports" directory created by DRAGEN.
    Accepts a Pooling file as argument, but use parent folder of `reports_dir`
    to look for the default 'SamplePools.csv' file.
    """
    args = parse_args()
    configure_logging(args.level)

    reports_dir = os.path.abspath(args.reports_dir)
    if args.pools is not None:
        pools_csv = os.path.abspath(args.pools)
    else:
        # Pooling file not provided, assuming that default file can be found 
        # in the parent folder of "Reports/", Dragen's default output folder
        #
        parent_dir = os.path.split(reports_dir)[0]
        pools_csv = f"{parent_dir}/SamplePools.csv"

    logging.info(f"Assembling XLSX Report from CSV files in output directory {reports_dir}")
    logging.debug(f"OPTIONS={args}")

    logging.info(f"Done.\n")


if __name__ == '__main__':
    sys.exit(main())

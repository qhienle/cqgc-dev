#!/usr/bin/env python3
"""
List samples in a given Run.

USAGE: list_run_samples.py RUN
       list_run_samples.py 20240130_LH00336_0009_A22GNV2LT3
       list_run_samples.py --help

Requires full name of the run as command-line argument.
Outputs a list of samples and families in file "samples_list.csv".
"""

import os, sys
import argparse
import logging
import json
import pandas as pd

# Set source path to CQGC-utils so that we can use relative imports
#
src_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(src_path)
from lib.nanuq import Nanuq


__version__ = 0.1


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="List samples in a given RUN.")
    parser.add_argument('run', help="Run ID for flowcell, ex: '20240130_LH00336_0009_A22GNV2LT3'")
    parser.add_argument('--logging-level', '-l', dest='level', default='info',
                        help="Logging level, can be 'debug', 'info', 'warning'. Default='info' [str]")
    return parser.parse_args()


def configure_logging(level):
    """
    Set logging level, based on the level names of the `logging` module.
    - level: [str] 'debug', 'info' or 'warning'
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


def get_nanuq_sample_data(cqgc_id):
    """
    Get from Nanuq family information for biosample `cqgc_id`.
    - `cqgc_id`: [str] sample identifier
    - Return: [dict]
    """
    sample_infos = {}
    nq = Nanuq()
    try:
        data = json.loads(nq.get_sample(cqgc_id))
    except Exception as e:
        logging.error(f"JSONDecodeError {e} could not decode biosample {cqgc_id}")
    else:
        logging.info(f"Got information for biosample {cqgc_id}")

    if len(data) == 0:
        logging.error(f"No samples retrieved from Nanuq for {cqgc_id}!")
    elif len(data) > 1:
        logging.debug(f"Number of samples retrieved from Nanuq is not 1.\n{data}")
    else: # len(data) == 1
        pass

    try:
        data[0]["patient"]["mrn"]
    except Exception as err:
        logging.warning(f"Could not find MRN for patient {cqgc_id}: {err}")
        data[0]["patient"]["mrn"] = '0000000'
    else:
        pass
    finally:
        sample_infos = {
            'sample_name': data[0]["ldmSampleId"],
            'biosample'  : data[0]["labAliquotId"],
            'relation'   : data[0]["patient"]["familyMember"],
            'gender'     : data[0]["patient"]["sex"],
            'ep_label'   : data[0]["patient"]["ep"],
            'mrn'        : data[0]["patient"]["mrn"],
            'status'     : data[0]["patient"]["status"],
            'family_id'  : data[0]["patient"].get("familyId", "-"),
            'birthdate'  : data[0]["patient"]["birthDate"],
            'project'    : data[0]["projectName"]
        }
    return sample_infos


def main(args):
    """
    For each sample listed in SampleSheet, collect metrics and generate reports
    - `args` : Command-line arguments, from `argparse`.
    - Returns: Lit of samples for Run in file 'samples_list.csv'
    """
    # Setup environment for this run. Results are written to folder "work_dir",
    # some information collected here will be used for case creation later on
    # Emedgene.
    #
    nq = Nanuq()

    fc_parts = nq.parse_run_name(args.run)
    fc_date  = fc_parts[0]
    fc_instr = fc_parts[1]
    fc_short = fc_parts[4]
    print(f"# Logging run {fc_parts}")
    
    # List samples. Maybe more precise to use the SampleSheet but not resilient
    #
    logging.info(f'Creating "samples_list.csv"')
    biosamples = []
    for tuple in nq.list_samples(args.run):
        biosamples.append(tuple[0])
    total = len(biosamples)
    logging.debug(f"Found {total} samples")

    # Collect family information from Nanuq for biosample (to build Case)
    # Build Pandas DataFrames from collected data for easier manipulations
    # Save DataFrame for samples to samples_list.csv, for later use
    #
    samples_families = [] # [{sample: val, gender: val, relation: val,...}, {...},...]
    for count, biosample in enumerate(biosamples, start=1):
        logging.info(f"Collecting family information for {biosample}, {count}/{total}")
        try:
            samples_families.append(get_nanuq_sample_data(biosample))
        except Exception as e:
            logging.error(f"In `get_nanuq_sample_data({biosample})`: {e}")
            logging.warning(f"COULD NOT RETRIEVE INFO FOR biosample {biosample}`. SKIPPING...")

    df_samples_families = pd.DataFrame(samples_families)
    df_samples_families = df_samples_families.sort_values(by=['family_id', 'relation'], ascending=[True, False])
    # Fix dates out of bounds with pd.Timestamp.min (eg: 11/11/1111) with errors='coerce'.
    # TODO: Check that downstream processes will accept null DateTime, 'NaT'.
    df_samples_families['birthdate'] = pd.to_datetime(df_samples_families['birthdate'], format='mixed', errors='coerce') # format='%d/%m/%Y')
    df_samples_families['flowcell_date'] = pd.to_datetime(fc_date, format='%Y-%m-%d')
    df_samples_families['flowcell'] = args.run

    df_samples_families.to_csv('samples_list.csv', index=None)
    logging.info(f"Collected family information into file 'samples_list.csv'")


if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)

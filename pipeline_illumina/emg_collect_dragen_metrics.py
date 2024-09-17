#!/usr/bin/env python3
"""
Collect DRAGEN metrics for samples in a given Run to generate reports in HTML
and CSV format. Data is collected from the file "{sample}.metrics.json" output
by DragenGermline (version 4+).

USAGE: emg_collect_dragen_metrics.py RUN
       emg_collect_dragen_metrics.py 20240130_LH00336_0009_A22GNV2LT3
       emg_collect_dragen_metrics.py --help

Requires full name of the run as command-ilne argument.
       
Produces reports in .html and .csv file, in a folder named 
{instrument}_{run_number} on `spxp-app02://staging2/dragen/`. For example:
spxp-app02://staging2/dragen/LH00336_0009
"""

import os, sys
import argparse
import logging
import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

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
    parser = argparse.ArgumentParser(description="Collect DRAGEN metrics for PRAGMatIQ samples in a given Run.")
    parser.add_argument('run', help="Run ID for flowcell, ex: '20240130_LH00336_0009_A22GNV2LT3'")
    parser.add_argument('--data-dir', '-d', dest = 'data_dir', 
                        help="Path to DragenGermline output. Ex: /staging/hiseq_raw/LH00336/20240130_LH00336_0009_A22GNV2LT3/Analysis/1/Data/DragenGermline")
    parser.add_argument('--samples-list-only', '-s', action='store_true', dest = 'samples_list_only', 
                        help="Only create the file 'samples_list.csv'")
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
        logging.warning(f"JSONDecodeError {e} could not decode biosample {cqgc_id}")
    else:
        logging.info(f"Got information for biosample {cqgc_id}")
    if len(data) != 1:
        logging.debug(f"Number of samples retrieved from Nanuq is not 1.\n{data}")
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
            'birthdate'  : data[0]["patient"]["birthDate"]
        }
    return sample_infos


def write_html_report(df, fc_short):
    """
    Write HTML report from data in `df`.
    - `df`: Pandas DataFrame.
    - Returns: HTML file.
    """
    fig1 = go.Figure(
        data = [
            go.Bar(name='Mapped Reads %', x=df['sample_name'], y=df['mapped_reads_pct'], yaxis='y', offsetgroup=1),
            go.Bar(name='SNPs (pass) %',  x=df['sample_name'], y=df['variants_snps_pass_pct'], yaxis='y2', offsetgroup=2),
            go.Bar(name='CNVs (amp+del)', x=df['sample_name'], y=df['cnvs_number'], yaxis='y3', offsetgroup=3)
        ],
        layout={
            'yaxis':  {'title': 'Percentage of mapped reads'},
            'yaxis2': {'title': 'Percentage of SNPs (pass)',  'overlaying': 'y', 'side': 'right'},
            'yaxis3': {'title': 'Number of CNVs (amplifications + deletions)', 'overlaying': 'y', 'side': 'left'},
            'barmode': 'group'
        }
    )
    #fig1.update_layout(barmode='group')

    fig2 = px.violin(df, y="cnv_coverage_uniformity", x="ep_label",
                     title="Coverage uniformity of CNVs per site",
                     color="ep_label", box=True, points="all", hover_data=df.columns)

    fig3 = px.scatter(df, x="cnv_coverage_uniformity", y="cnvs_number", 
                      title="Number of CNVs vs Coverage uniformity",
                      hover_data=['sample_name', 'cnv_coverage_uniformity', 'ep_label'], color='ep_label')
    fig3.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    fig4 = px.scatter(df, x="uniformity_of_coverage_pct_gt_02mean_over_genome", y="variants_snps_pass_pct", 
                      title="Number of variants passing filter and uniformity of coverage >0.2 mean over genome",
                      hover_data=['sample_name', 'pct_of_genome_with_coverage_20x_inf', 'uniformity_of_coverage_pct_gt_02mean_over_genome', 'ep_label'], 
                      color='ep_label')
    fig4.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    fig5 = px.scatter(df, x="average_alignment_coverage_over_genome", y="cnv_coverage_uniformity", 
                      title="CNV coverage uniformity vs Average coverage",
                      hover_data=['sample_name', 'cnv_coverage_uniformity', 'average_alignment_coverage_over_genome','ep_label'], 
                      color='ep_label')             
    fig5.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    css = """
            html body {
                background: #f8f9f9;
                font-size: 11px;
            }
            table { 
                background: white;
                width: 100%;
                max-width: 800px;
                margin: 50px auto;
                border-collapse: collapse; 
                border-spacing: 1; 
                border-radius: 6px;
                overflow: hidden;
                position: relative;
                font-family: Arial;
                font-size: 11px;
            }
            /* Zebra striping */
            tr:nth-of-type(odd) { 
                background: #eee; 
            }
            th { 
                background: #3498db; 
                color: white; 
                font-weight: bold; 
            }
            td, th { 
                padding: 10px; 
                border: 1px solid #ccc; 
                text-align: left; 
                font-size: 18px;
            }
    """
    out_html = f"{fc_short}_metrics.html"
    title    = f"{fc_short} Samples Metrics"

    with open(out_html, 'w') as fh:
        fh.write('<!doctype html>\n<html>\n\t<head>\n')
        fh.write(f'\t\t<title>{title}</title>\n\t\t<meta charset="UTF-8">\n')
        fh.write(f'\t\t<style type="text/css">{css}\t\t</style>\n')
        fh.write('\t</head>\n\t<body>\n')
        fh.write(f'\t\t<h1>{fc_short}</h1>\n\t\t')
        fh.write(df.to_html(index=False))
        fh.write(fig1.to_html(full_html=False, include_plotlyjs='cdn'))
        fh.write(fig2.to_html(full_html=False, include_plotlyjs='cdn'))
        fh.write(fig3.to_html(full_html=False, include_plotlyjs='cdn'))
        fh.write(fig4.to_html(full_html=False, include_plotlyjs='cdn'))
        fh.write(fig5.to_html(full_html=False, include_plotlyjs='cdn'))
        fh.write('\n\t</body>\n</html>')


def main(args):
    """
    For each sample listed in SampleSheet, collect metrics and generate reports
    - `args` : Command-line arguments, from `argparse`.
    - Returns: Two files named ./archives_metrics.csv and archives_metrics.html
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
    work_dir = f"/staging2/dragen/{fc_short}"
    print(f"# Logging run {fc_parts}")
    
    try:
        os.mkdir(work_dir)
    except FileExistsError as err:
        logging.info(err)
    finally:
        os.chdir(work_dir)

    # List samples. Maybe more precise to use the SampleSheet's [DragenGermline]
    # section, but not very resilient.
    #
    logging.info(f'Creating "samples_list.csv"')
    biosamples = []
    for tuple in nq.list_samples(fc_short):
        biosamples.append(tuple[0])
    total = len(biosamples)
    logging.debug(f"Found {total} samples")

    # Build Pandas DataFrames from collected data for easier manipulations

    # Collect family information from Nanuq for biosample (to build Case)
    # Save DataFrame for samples to samples_list.csv, for later use
    #
    samples_families = [] # [{sample: val, gender: val, relation: val,...}, {...},...]
    for count, biosample in enumerate(biosamples, start=1):
        logging.info(f"Collecting family information for {biosample}, {count}/{total}")
        samples_families.append(get_nanuq_sample_data(biosample))

    df_samples_families = pd.DataFrame(samples_families)
    df_samples_families = df_samples_families.sort_values(by=['family_id', 'relation'], ascending=[True, False])
    df_samples_families['birthdate'] = pd.to_datetime(df_samples_families['birthdate'], format='mixed') # format='%d/%m/%Y')
    df_samples_families['flowcell_date'] = pd.to_datetime(fc_date, format='%Y%m%d')
    df_samples_families['flowcell'] = args.run
    df_samples_families.to_csv('samples_list.csv', index=None)
    logging.info(f"Collected family information into file 'samples_list.csv'")

    sys.exit() if args.samples_list_only else None

    # Collect samples metrics from DragenGermline analyses
    # Save DataFrame so that we can merge metrics with family information
    #
    samples_metrics  = {} # {biosample: {metric1: value, metric2: value2, ...}}

    if args.data_dir is not None:
        data_dir = args.data_dir
    else:
        # Check if there is more than one dir under "./Analysis/". If yes, 
        # we ask to use --data-dir so that we don't have to guess.
        #
        analyses = os.listdir()
        if len(analyses) == 1 and analyses[0] == 1:
            data_dir = f"/staging/hiseq_raw/{fc_instr}/{args.run}/Analysis/1/Data/DragenGermline"
        else:
            logging.info(f"More than one Analysis directory found.")
            sys.exit('Please specify which one to use with option `--data-dir`')

    for count, biosample in enumerate(biosamples, start=1):
        logging.info(f"Collecting family information for {biosample}, {count}/{total}")
        metrics_file = f"{data_dir}/{biosample}/germline_seq/{biosample}.metrics.json"
        try:
            with open(metrics_file, 'r') as fh:
                metrics = json.load(fh)['Attributes']['illumina_dragen_complete_v0_4']
        except Exception as e:
            logging.warn(f"While opening file {metrics_file} for {biosample}, got ERROR: {e}")
            samples_metrics[biosample] = None
        else:
            samples_metrics[biosample] = metrics
            logging.info(f"Collected metrics for {biosample}")

    df_samples_metrics  = pd.DataFrame.from_dict(samples_metrics, orient="index")
    df_samples_metrics['biosample'] = df_samples_metrics.index
    df_samples_metrics['cnvs_number'] = df_samples_metrics['cnv_number_of_amplifications'] + df_samples_metrics['cnv_number_of_deletions']
    logging.info(f"Collected DragenGermline metrics")

    # Subset columns for report, join infos for sample name and site label
    #
    subset_cols =  ['biosample',
                   'mapped_reads_pct', 
                   'average_alignment_coverage_over_genome',
                   'variants_snps_pass_pct',
                   'cnvs_number',
                   'cnv_number_of_amplifications', 
                   'cnv_number_of_passing_amplifications',
                   'cnv_number_of_deletions',
                   'cnv_number_of_passing_deletions',
                   'cnv_coverage_uniformity',
                   'pct_of_genome_with_coverage_20x_inf',
                   'uniformity_of_coverage_pct_gt_02mean_over_genome',
                   'mean_median_autosomal_coverage_ratio_over_genome',
                   'number_unique_mapped_reads_excl_duplicates',
                   'number_unique_mapped_reads_excl_duplicates_pct',
                   'number_of_duplicate_marked_reads'
                   ]
    df_subset_metrics = df_samples_metrics[subset_cols]
    df_report = pd.merge(df_subset_metrics, df_samples_families[['biosample', 'sample_name', 'ep_label']], on='biosample', how="outer")
    logging.info(f"Writing report for {fc_short}:\n{df_report}")
    df_report.to_csv(f'{fc_short}_metrics.csv', index=None)
    write_html_report(df_report, fc_short)


if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)

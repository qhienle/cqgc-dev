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


def list_dragengermline_samples(samplesheet):
    """
    Get list of samples from [DragenGermline] section in Nanuq a SampleSheet
    - `samplesheet`: [str] path to file SampleSheet.csv.
    - Returns: [list] of samples
    """
    samples = []
    try:
        with open(samplesheet, 'r') as fh:
            content = fh.read().splitlines()
    except FileNotFoundError as err:
        logging.error(f"Could not open {samplesheet}: {err}")
    else:
        for line in content:
            if line.startswith('['):
                section = line
            else:
                if section.startswith('[DragenGermline_Data]'):
                    cols = line.split(',')
                    if not line.startswith('Sample_ID') and len(cols) > 1:
                        samples.append(cols[0])
    return samples


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
            'label'      : data[0]["patient"]["ep"],
            'mrn'        : data[0]["patient"]["mrn"],
            'cohort_type': data[0]["patient"]["designFamily"],
            'status'     : data[0]["patient"]["status"],
            'Family Id'  : data[0]["patient"].get("familyId", "-"),
            'date_of_birth(YYYY-MM-DD)': data[0]["patient"]["birthDate"]
        }
        #TODO: Add 'pid', 'phenotypes', 'hpos', 'filenames'
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
            'yaxis3': {'title': 'Number of CNVs (aplifications + deletions)', 'overlaying': 'y', 'side': 'left'},
            'barmode': 'group'
        }
    )
    #fig1.update_layout(barmode='group')

    fig2 = px.violin(df, y="cnv_coverage_uniformity", x="label",
                     title="Coverage uniformity of CNVs per site",
                     color="label", box=True, points="all", hover_data=df.columns)

    fig3 = px.scatter(df, x="cnv_coverage_uniformity", y="cnvs_number", 
                      title="Number of CNVs vs Coverage uniformity",
                      hover_data=['sample_name', 'cnv_coverage_uniformity', 'label'], color='label')
    fig3.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    fig4 = px.scatter(df, x="uniformity_of_coverage_pct_gt_02mean_over_genome", y="variants_snps_pass_pct", 
                      title="Number of variants passing filter and uniformity of coverage >0.2 mean over genome",
                      hover_data=['sample_name', 'pct_of_genome_with_coverage_20x_inf', 'uniformity_of_coverage_pct_gt_02mean_over_genome', 'label'], 
                      color='label')
    fig4.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    fig5 = px.scatter(df, x="average_alignment_coverage_over_genome", y="cnv_coverage_uniformity", 
                      title="CNV coverage uniformity vs Average coverage",
                      hover_data=['sample_name', 'cnv_coverage_uniformity', 'average_alignment_coverage_over_genome','label'], 
                      color='label')             
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
    nq   = Nanuq()

    fc_parts = nq.parse_run_name(args.run)
    fc_date  = fc_parts[0]
    fc_instr = fc_parts[1]
    fc_short = fc_parts[4]
    print(f"# Logging run {fc_parts}")
    
    if args.data_dir is not None:
        data_dir = args.data_dir
    else:
        # If there is more than one dir under "./Analysis/", use --data-dir"
        data_dir = f"/staging/hiseq_raw/{fc_instr}/{args.run}/Analysis/1/Data/DragenGermline"
    work_dir = f"/staging2/dragen/{fc_short}"
    try:
        os.mkdir(work_dir)
    except FileExistsError as err:
        logging.info(err)
    finally:
        os.chdir(work_dir)

    # Collect data from "*.metrics.json" for all the samples into a Pandas 
    # DataFrame. List of samples is taken from [DragenGermline_Data] section of
    # "SampleSheet.csv". NOTE: SampleSheet from Nanuq GET API doesn't yet have 
    # [DragenGermline_Data] section.
    #
    biosamples = list_dragengermline_samples(f"{data_dir}/SampleSheet.csv")
    total      = len(biosamples)

    samples_metrics  = {} # {biosample: {metric1: value, metric2: value2, ...}}
    samples_families = [] # [{sample: val, gender: val, relation: val,...}, {...},...]

    for count, biosample in enumerate(biosamples, start=1):
        logging.info(f"Processing {biosample}, {count}/{total}")

        # Collect metrics for this biosample
        #
        metrics_file = f"{data_dir}/{biosample}/germline_seq/{biosample}.metrics.json"
        with open(metrics_file, 'r') as fh:
            metrics = json.load(fh)['Attributes']['illumina_dragen_complete_v0_4']
        samples_metrics[biosample] = metrics
        logging.info(f"Collected metrics for {biosample}")

        # Collect family information from Nanuq for biosample (to build Case)
        #
        samples_families.append(get_nanuq_sample_data(biosample))

    # Build Pandas DataFrames from collected data for easier manipulations
    #
    df_samples_families = pd.DataFrame(samples_families)
    df_samples_families = df_samples_families.sort_values(by=['Family Id', 'relation'], ascending=[True, False])

    df_samples_metrics  = pd.DataFrame.from_dict(samples_metrics, orient="index")
    df_samples_metrics['biosample'] = df_samples_metrics.index
    df_samples_metrics['cnvs_number'] = df_samples_metrics['cnv_number_of_amplifications'] + df_samples_metrics['cnv_number_of_deletions']
    logging.info(f"Built DataFrames from metrics and family information")

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
    df_report = pd.merge(df_subset_metrics, df_samples_families[['biosample', 'sample_name', 'label']], on='biosample', how="outer")
    logging.info(f"Writing report for {fc_short}:\n{df_report}")
    write_html_report(df_report, fc_short)


if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)

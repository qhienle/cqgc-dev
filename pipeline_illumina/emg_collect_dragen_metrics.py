#!/usr/bin/env python3
"""
Collect DRAGEN metrics for PRAGMatIQ samples in a given Run.

USAGE: emg_collect_dragen_metrics.py RUN
       emg_collect_dragen_metrics.py --help

Lorem ipsum
"""

import os
import argparse
import logging
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from glob import glob

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


def list_dragengermline_samples(samplesheet='SampleSheet.csv'):
    """
    Get list of samples from [DragenGermline] section in Nanuq a SampleSheet
    - `samplesheet`: [str] path to file SampleSheet.csv.
    - Returns: [list] of samples
    """
    samples = []
    with open(samplesheet, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('['):
                section = line
            else:
                if section.startswith('[DragenGermline_Data]'):
                    cols = line.split(',')
                    if not line.startswith('Sample_ID') and len(cols) > 1:
                        samples.append(cols[0])
    return samples


def write_html_report(df, fc_short):
    """
    Write HTML report from data in `df`.
    - `df`: Pandas DataFrame.
    - Returns: HTML file.
    """
    fig1 = go.Figure(
        data = [
            go.Bar(name='NumOfCNVs',  x=df['Sample'], y=df['NumOfCNVs'], yaxis='y', offsetgroup=1),
            go.Bar(name='NumOfSNPs',  x=df['Sample'], y=df['NumOfSNPs'], yaxis='y2', offsetgroup=2),
            go.Bar(name='NumOfReads', x=df['Sample'], y=df['NumOfReads'], yaxis='y3', offsetgroup=3)
        ],
        layout={
            'yaxis': {'title': 'NumOfCNVs'},
            'yaxis2': {'title': 'NumOfSNPs',  'overlaying': 'y', 'side': 'right'},
            'yaxis3': {'title': 'NumOfReads', 'overlaying': 'y', 'side': 'left'},
            'barmode': 'group'
        }
    )
    #fig1.update_layout(barmode='group')

    fig2 = px.violin(df, y="Coverage uniformity", x="Site",
                     title="Coverage uniformity of CNVs per site",
                     color="Site", box=True, points="all", hover_data=df.columns)

    fig3 = px.scatter(df, x="Coverage uniformity", y="NumOfCNVs", 
                      title="Number of CNVs vs Coverage uniformity",
                      hover_data=['Sample', 'Coverage uniformity', 'Site'], color='Site')
    fig3.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    fig4 = px.scatter(df, x="Uniformity of coverage (PCT > 0.4*mean) over genome", y="NumOfSNPs", 
                      title="Number of SNPs and Uniformity of coverage",
                      hover_data=['Sample', 'PCT coverage >20x', 'Uniformity of coverage (PCT > 0.4*mean) over genome', 'Site'], 
                      color='Site')
    fig4.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=2, color='DarkSlateGrey')), selector=dict(mode='markers'))

    fig5 = px.scatter(df, x="Average coverage", y="Coverage uniformity", 
                      title="Coverage uniformity vs Average coverage",
                      hover_data=['Sample', 'Coverage uniformity', 'Average coverage','Site'], 
                      color='Site')             
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
    Get list of samples from the SampleSheet.
    For each sample, collect metrics and generate a report
    - `args` : Command-line arguments, from `argparse`.1
    - Returns: A CSV file named `./archives_metrics.csv`.
    """
    # Setup environment for this run
    #
    date, instr, run_nb, flowcell = args.run.split('_')
    fc_short = f"{instr}_{run_nb}"
    work_dir = f"/staging2/dragen/{fc_short}"
    if args.data_dir is not None:
        data_dir = args.data_dir
    else:
        data_dir = f"/staging/hiseq_raw/{instr}/{args.run}/Analysis/1/Data/DragenGermline"
    os.mkdir(work_dir)
    os.chdir(work_dir)

    # Process list of samples contained in samples_list file.
    #
    samples = list_dragengermline_samples()
    total   = len(samples)

    # Initialize empty Pandas DataFrames
    #
    df_metrics   = get_metrics_from_log(None)
    df_coverages = get_coverage_metrics(None)
    
    for count, sample in enumerate(samples, start=1):
        logging.info(f"Processing {sample}, {count}/{total}")

        df_metrics = pd.concat([df_metrics, get_metrics_from_log(sample, log_files=log_files)], ignore_index=True)
        df_coverages   = pd.concat([df_coverages, get_coverage_metrics(sample, coverage_files=coverage_files)], ignore_index=True)

    df_metrics.drop_duplicates(inplace=True)
    df_metrics.reset_index(inplace=True)
    df_coverages.drop_duplicates(inplace=True)
    df_coverages.reset_index(inplace=True)
    df = df_metrics.merge(df_coverages, on='Sample', how='outer')
    df.drop(['index_x'], axis=1, inplace=True)

    df1 = df[['Sample', 'NumOfReads', 'NumOfSNPs', 'NumOfCNVs',
        'Average coverage', 'Coverage uniformity',
        'CNV Number of amplifications', 'CNV Number of passing amplifications', 
        'CNV Number of deletions', 'CNV Number of passing deletions', 
        'PCT coverage >20x',
        'Uniformity of coverage (PCT > 0.2*mean) over genome', 
        'Uniformity of coverage (PCT > 0.4*mean) over genome', 
        'Mean/Median autosomal coverage ratio over genome',
        'Number of unique & mapped reads',
        'Number of unique & mapped reads PCT',
        'Number of duplicate marked reads PCT', 
        'Percent Autosome Callability']] #, 'Estimated sample contamination']]
    df1.drop_duplicates(inplace=True)
    df1 = df1.astype({'NumOfReads': 'Int64', 'NumOfSNPs': 'Int64', 'NumOfCNVs': 'Int64',
                    'Average coverage': 'Float64', 'Coverage uniformity': 'Float64',
                    'CNV Number of amplifications': 'Int64', 'CNV Number of passing amplifications': 'Int64', 
                    'CNV Number of deletions': 'Int64', 'CNV Number of passing deletions': 'Int64', 
                    'PCT coverage >20x': 'Float64',
                    'Uniformity of coverage (PCT > 0.2*mean) over genome': 'Float64', 
                    'Uniformity of coverage (PCT > 0.4*mean) over genome': 'Float64', 
                    'Mean/Median autosomal coverage ratio over genome': 'Float64',
                    'Number of unique & mapped reads': 'Int64',
                    'Number of unique & mapped reads PCT': 'Float64',
                    'Number of duplicate marked reads PCT': 'Float64', 
                    'Percent Autosome Callability': 'Float64'})
    logging.info(f"Current Metrics DataFrame:\n{df1}")

    # Combine dataframe with samples_list.csv and generate figures for the HTML report
    #
    df_samples = pd.read_csv('samples_list.csv', encoding="latin-1")
    df = df1.merge(df_samples, how='inner')
    df.to_csv(f'{args.run}_metrics.csv', index=None)
    write_html_report(df, args.run)


if __name__ == '__main__':
    args = parse_args()
    configure_logging(args.level)
    main(args)
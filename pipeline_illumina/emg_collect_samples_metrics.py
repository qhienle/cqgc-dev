#!/usr/bin/env python3
"""
Collect metrics for PRAGMatIQ samples.

USAGE: emg_get_samples_metrics.py SAMPLES
       emg_get_samples_metrics.py --help

Parse Emedgene's logs files to get analysis metrics for PRAGMatIQ samples 
listed in file SAMPLES, a CSV file in which the names of samples are in the 1st
column. Ex: `samples_list.csv` output from `emg_make_batch_from_nanuq.py`:

    Sample,CQGC_ID,Site,Date
    GM231651,22293,CHUSJ,2023-08-09
    23-05982-T1,22282,CHUS,2023-08-09
    3042652455,22256,CHUQ,2023-08-09
    (...)

-d|--directory: Folder where log files are chached. If the folder is already
present, log file therein are re-used, instead of downloading from AWS S3.
Default folder is `emg_logs`.
    
Metrics collected are written to a CSV table and an HTML report file. 
Log files are downloaded from EMG's S3 bucket using the CLI `aws` in a tmp dir.  
"""

import os
import argparse
import logging
import subprocess
import pandas as pd
from glob import glob

__version__ = 0.1


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Collect metrics for PRAGMatIQ samples")
    parser.add_argument('-s', '--samples', default="samples_list.csv", 
                        help="Filename to CSV list of samples. Default=`samples_list.csv` [str]")
    parser.add_argument('-d', '--directory', dest='dir', default="emg_logs", 
                        help="Directory containing EMG log files. Default='emg_logs' [str]")
    parser.add_argument('-l', '--logging-level', dest='level', default='info',
                        help="Logging level, can be 'debug', 'info', 'warning'. Default='info' [str]")
    parser.add_argument('-p', '--profile', default='emedgene',
                        help="Emedgene profile, can be 'emedgene' or 'emedgene-eval'. Default='emedgene' [str]")
    return(parser.parse_args())


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


def get_samples_list(file):
    """
    Get list of samples to archive from a CSV `file`.
    - `file`: [str] path to CSV file with list of sample names in 1st column.
    - Returns: [list] of samples
    """
    samples = []
    with open(file, 'r') as fh:
        next(fh) # Skip header
        for line in fh.readlines():
            sample = line.split(',')[0]
            samples.append(sample)
    return samples


def download_emg_s3_logs(sample, profile='emedgene', logsdir='emg_logs'):
    """
    List files available on AWS bucket for Emedgene `profile` and download 
    to logsdir three files for collecting dragen metrics of `sample`: 
    '_sample.log', '.dragen.bed_coverage_metrics.csv' and '.dragen.cnv.vcf.gz'
    - sample : Name of sample to retrive [str]
    - profile: AWS credentials, either `emedgene` or `emedgene-eval` [str]
    - returns: Downloaded files under folder `logsdir` [list]
        NOTE: As there may be more than one of each "metrics" and "vcf" files.
        these are tagged with the version hash of the subdir on s3.
    """

    os.mkdir(logsdir) if not os.path.isdir(logsdir) else None
    site = 's3://cac1-prodca-emg-auto-results'
    domains = {'emedgene': 'CHU_Sainte_Justine', 'emedgene-eval': 'Ste_Justine_eval'}
    url = f"{site}/{domains[profile]}/{sample}"
        
    ls = subprocess.run(['aws', 's3', '--profile', profile, 'ls', '--recursive', url], capture_output=True, text=True)
    files = []
    for line in ls.stdout.splitlines():
        s3_file    = line.split()[-1]
        s3_url     = f"{site}/{s3_file}"
        path_parts = s3_file.split('/')
        file       = path_parts[-1]
        if file.endswith("_sample.log"):
            if not os.path.isfile(f"{logsdir}{os.sep}{file}"):
                logging.info(f"Log file {logsdir}{os.sep}{file} not found. Downloading from S3...")
                subprocess.run(['aws', 's3', '--profile', profile, 'cp', s3_url, logsdir], check=True)
            files.append(f"{logsdir}{os.sep}{file}")
        elif file.endswith(".dragen.bed_coverage_metrics.csv"):
            output = f"{logsdir}{os.sep}{sample}_{path_parts[4]}.dragen.bed_coverage_metrics.csv"
            if not os.path.isfile(output):
                logging.info(f"Log file {output} not found. Downloading from S3...")
                subprocess.run(['aws', 's3', '--profile', profile, 'cp', s3_url, output], check=True)
            files.append(output)
        elif file.endswith(".dragen.cnv.vcf.gz"):
            output = f"{logsdir}{os.sep}{sample}_{path_parts[4]}.dragen.cnv.vcf.gz"
            if not os.path.isfile(output):
                logging.info(f"Log file {output} not found. Downloading from S3...")
                subprocess.run(['aws', 's3', '--profile', profile, 'cp', s3_url, output], check=True)
            files.append(output)
    return files


def glob_files(pattern):
    """
    Glob list of files from pattern and checks that list is not empty
    """
    files = glob(pattern)
    if len(files) == 0:
        logging.debug(f"More than one log file found with pattern {pattern}")
    elif len(files) == 0:
        logging.warning(f"No file found with pattern {pattern}")
    return files


def get_metrics_from_log(sample, log_files=[]):
    """
    Get metrics from log file. Metrics for "Percent of genome coverage over 20x" 
    and "Average coverage" are not consistently formatted. Get these infos from
    file *.dragen.bed_coverage_metrics.csv
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame, with the following information per sample
        - Number of reads
        - SNPs
        - Average coverage for CNV
        - CNV CoverageUniformity
        - Percent Autosome Callability
    """
    metrics = [] # [[Sample, Log filename, Number of reads, SNPs, CNV Average coverage, Coverage uniformity], [],...]
    logging.debug(f"List of logfiles to parse: {log_files}")
    for log in log_files:
        logname = os.path.basename(log)
        with open(log, "r") as fh:
            lines = fh.readlines()
        reads = ''
        snps  = ''
        amplifications      = ''
        pass_amplifications = ''
        deletions           = '' 
        pass_deletions      = ''
        coverage_uniformity = ''
        callability         = ''
        contamination       = '' 
        mapped_reads        = '' 
        mapped_reads_pct    = '' 
        duplicate_reads_pct = '' 
        for line in lines:
            line_parts = line.rstrip().split()

            # For some reason, there are often duplicate lines for these metrics in the log file, 
            # but not in the original bed_coverage_metrics.csv file.
            # The information is sometimes presented with the bytestring (b'blah...\n') symbols 
            # and need to be reformatted.
            #
            if 'Number of reads:' in line:
                reads = line_parts[-1].replace("\\n'", "")
            # CNV?
            # elif 'Average alignment coverage over genome' in line and 'COVERAGE SUMMARY' in line:
            #     cnv_avg_coverage = line_parts[-1].replace("\\n'", "")
            elif 'Coverage uniformity' in line:
                coverage_uniformity = line_parts[-1].replace("\\n'", "")
            elif 'Number of amplifications' in line and 'CNV SUMMARY' in line:
                amplifications = int(line_parts[-1].replace("\\n'", ""))
            elif 'Number of passing amplifications' in line and 'CNV SUMMARY' in line:
                pass_amplifications = int(line_parts[-2].replace("\\n'", ""))
            elif 'Number of deletions' in line and 'CNV SUMMARY' in line:
                deletions = int(line_parts[-1].replace("\\n'", ""))
            elif 'Number of passing deletions' in line and 'CNV SUMMARY' in line:
                pass_deletions = int(line_parts[-2].replace("\\n'", ""))
            elif 'SNPs' in line and line_parts[11] == "SNPs":
                snps = line_parts[12].replace("\\n'", "")
            elif 'Percent Autosome Callability' in line:
                callability = line_parts[14].replace("\\n'", "")
            elif 'Number of unique & mapped reads' in line and 'MAPPING/ALIGNING SUMMARY' in line:
                mapped_reads     = line_parts[-2]
                mapped_reads_pct = line_parts[-1].replace("\\n'", "")
            elif 'Number of duplicate marked reads' in line and 'MAPPING/ALIGNING SUMMARY' in line:
                duplicate_reads_pct = line_parts[-1].replace("\\n'", "")
            elif 'Estimated sample contamination' in line:
                if line_parts[12] != 'standard':
                    contamination = line_parts[12].replace("\\n'", "")
        cnvs = amplifications + deletions
        metrics.append([sample, logname, reads, snps, cnvs, coverage_uniformity, 
                        amplifications, pass_amplifications, deletions, pass_deletions, 
                        mapped_reads, mapped_reads_pct, duplicate_reads_pct,
                        callability, contamination])
    df = pd.DataFrame(metrics, columns = ["Sample", "Log filename", 
                                            "NumOfReads", "NumOfSNPs", "NumOfCNVs",
                                            "Coverage uniformity", 
                                            "CNV Number of amplifications", 
                                            "CNV Number of passing amplifications", 
                                            "CNV Number of deletions", 
                                            "CNV Number of passing deletions", 
                                            'Number of unique & mapped reads',
                                            'Number of unique & mapped reads PCT',
                                            'Number of duplicate marked reads PCT', 
                                            "Percent Autosome Callability", 
                                            "Estimated sample contamination"])
    return df


def get_coverage_metrics(sample, coverage_files=[]):
    """
    Get coverage metrics for `sample`.
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame, with the following information per sample
        - average coverage
        - PCT coverage >20x
        - Uniformity of coverage (PCT > 0.2*mean) over genome
    """
    coverages = []
    logging.debug(f"List of logfiles to parse: {coverage_files}")

    for file in coverage_files:
        path_parts   = os.path.split(file)
        version      = os.path.basename(path_parts[0])
        avg_coverage = ''
        coverage_20x = ''
        uniformity_coverage_02 = ''
        uniformity_coverage_04 = ''
        autosomal_cover_ratio  = ''
        with open(file, "r") as fh:
            for line in fh:
                cols = line.rstrip().split(',')
                if cols[2].startswith('Average alignment coverage over genome'):
                    avg_coverage = cols[3]
                elif '20x: inf' in cols[2]:
                    coverage_20x = cols[3]
                    # NB: Different versions of DRAGEN can have more or less whitespaces
                    # between the left bracket '[' and '20x'. For example:
                    # v1.2.2_dragen3.9.5-hg38_2bd1884/GM230658.dragen.bed_coverage_metrics.csv: 
                    # COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: inf),92.17
                    # v1.2.2_dragen4.0.3-hg38_0edf29a/GM230732.dragen.bed_coverage_metrics.csv:
                    # COVERAGE SUMMARY,,PCT of genome with coverage [  20x: inf),80.13
                elif 'Uniformity of coverage' in cols[2]:
                    if '(PCT > 0.2*mean)' in cols[2]:
                        uniformity_coverage_02 = cols[3]
                    elif '(PCT > 0.4*mean)' in cols[2]:
                        uniformity_coverage_04 = cols[3]
                elif 'Mean/Median autosomal coverage ratio over genome' in cols[2]:
                        autosomal_cover_ratio = cols[3]
        coverages.append([sample, version, avg_coverage, coverage_20x, uniformity_coverage_02, uniformity_coverage_04, autosomal_cover_ratio])
    return pd.DataFrame(coverages, columns=['Sample', 
                                            'Log coverage', 
                                            'Average coverage', 
                                            'PCT coverage >20x', 
                                            'Uniformity of coverage (PCT > 0.2*mean) over genome', 
                                            'Uniformity of coverage (PCT > 0.4*mean) over genome',
                                            'Mean/Median autosomal coverage ratio over genome'])


def main(args):
    """
    From a list of samples, retrieve several metrics.
    - `args` : Command-line arguments, from `argparse`.1
    - Returns: A CSV file named `./archives_metrics.csv`.
    """
    configure_logging(args.level)
    workdir = os.getcwd()
    logsdir = args.dir
    if os.path.isdir(logsdir):
        logging.info(f"Logs directory, '{logsdir}', already exists")
    else:
        try: 
            os.mkdir(logsdir)
        except FileNotFoundError as e:
            logging.error(f"{e}; logsdir={logsdir}")
    
    # Initialize empty Pandas DataFrames
    #
    df_metrics   = get_metrics_from_log(None)
    df_coverages = get_coverage_metrics(None)
    
    # Process list of samples contained in samples_list file.
    #
    samples = get_samples_list(args.samples)
    total   = len(samples)
    for count, sample in enumerate(samples, start=1):
        logging.info(f"Processing {sample}, {count}/{total}")
        logs = download_emg_s3_logs(sample, profile=args.profile, logsdir=logsdir)
        logging.debug(f"S3 logfiles: {logs}")

        log_files  = glob_files(f"{args.dir}/{sample}_v*_sample.log")
        df_metrics = pd.concat([df_metrics, get_metrics_from_log(sample, log_files=log_files)], ignore_index=True)

        coverage_files = glob_files(f"{args.dir}/{sample}*.dragen.bed_coverage_metrics.csv")
        df_coverages   = pd.concat([df_coverages, get_coverage_metrics(sample, coverage_files=coverage_files)], ignore_index=True)

    df_metrics.drop_duplicates(inplace=True)
    df_coverages.drop_duplicates(inplace=True)

    df_metrics.reset_index(inplace=True)
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
    print(df1[['NumOfCNVs', 'CNV Number of amplifications', 'CNV Number of deletions']])

    df1_csv = workdir + os.sep + 'archives_metrics.csv'
    df1.to_csv(df1_csv, index=None)
    logging.info(f"Current Metrics DataFrame:\n{df1}")


def _test(args):
    samples = get_samples_list(args.samples)
    print(samples)


if __name__ == '__main__':
    args = parse_args()
    main(args)
    #_test(args)

#!/usr/bin/env python3
"""
Collect metrics for PRAGMatIQ samples.

USAGE: emg_get_samples_metrics.py SAMPLES
       emg_get_samples_metrics.py --help

Parse Emedgene's logs files to get analysis metrics for PRAGMatIQ samples 
listed in file SAMPLES, a CSV file in which the names of samples are in the 1st
column. Ex: `samples_list.txt` output from `emg_make_batch_from_nanuq.py`:

    sample_name,biosample,label,fc_date
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
import gzip
import pandas as pd
from glob import glob

__version__ = 0.1


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Collect metrics for PRAGMatIQ samples")
    parser.add_argument('-s', '--samples', default="samples_list.txt", 
                        help="Filename to CSV list of samples. Default=`samples_list.txt` [str]")
    parser.add_argument('-d', '--directory', dest='dir', default="emg_logs", 
                        help="Directory containing EMG log files. Default='emg_logs' [str]")
    parser.add_argument('-l', '--logging-level', dest='level', default='info',
                        help="Logging level, can be 'debug', 'info', 'warning'. Default='info' [str]")
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
    selected log files required to collect metrics for `sample`.
    - sample:  Name of sample to retrives [str]
    - profile: AWS credentials, either `emedgene` or `emedgene-eval` [str]
    - returns: Downloaded files under folder `./emg_logs/` [list]
    """
    os.mkdir(logsdir) if not os.path.isdir(logsdir) else None
    site = 's3://cac1-prodca-emg-auto-results'
    domains = {'emedgene': 'CHU_Sainte_Justine', 'emedgene-eval': 'Ste_Justine_eval'}
    url = f"{site}/{domains[profile]}/{sample}"
        
    ls = subprocess.run(['aws', 's3', '--profile', profile, 'ls', '--recursive', url], capture_output=True, text=True)
    files = []
    for line in ls.stdout.splitlines():
        file    = line.split()[-1].split('/')[-1]
        s3_file = f"{site}/{line.split()[-1]}"
        if file.endswith("_sample.log"):
            if not os.path.isfile(f"{logsdir}{os.sep}{file}"):
                print(f"Log file {logsdir}{os.sep}{file} not found. Downloading from S3...")
                subprocess.run(['aws', 's3', '--profile', profile, 'cp', s3_file, logsdir], check=True)
            files.append(file)
        elif line.endswith(".dragen.bed_coverage_metrics.csv"):
            if not os.path.isfile(f"{logsdir}{os.sep}{file}"):
                print(f"Log file {logsdir}{os.sep}{file} not found. Downloading from S3...")
                subprocess.run(['aws', 's3', '--profile', profile, 'cp', s3_file, logsdir], check=True)
            files.append(file)
        elif line.endswith(".dragen.cnv.vcf.gz"):
            if not os.path.isfile(f"{logsdir}{os.sep}{file}"):
                print(f"Log file {logsdir}{os.sep}{file} not found. Downloading from S3...")
                subprocess.run(['aws', 's3', '--profile', profile, 'cp', s3_file, logsdir], check=True)
            files.append(file)
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


def get_metrics_from_log(sample, logfiles=[]):
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
    logging.debug(f"List of logfiles to parse: {logfiles}")
    for log in logfiles:
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
                amplifications = line_parts[-1].replace("\\n'", "")
            elif 'Number of passing amplifications' in line and 'CNV SUMMARY' in line:
                pass_amplifications = line_parts[-2].replace("\\n'", "")
            elif 'Number of deletions' in line and 'CNV SUMMARY' in line:
                deletions = line_parts[-1].replace("\\n'", "")
            elif 'Number of passing deletions' in line and 'CNV SUMMARY' in line:
                pass_deletions = line_parts[-2].replace("\\n'", "")
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
        metrics.append([sample, logname, reads, snps, coverage_uniformity, 
                        amplifications, pass_amplifications, deletions, pass_deletions, 
                        mapped_reads, mapped_reads_pct, duplicate_reads_pct,
                        callability, contamination])
    return pd.DataFrame(metrics, columns = ["Sample", "Log filename", 
                                            "NumOfReads", "NumOfSNPs",
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


def get_coverage_metrics(sample):
    """
    Get coverage metrics for `sample`.
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame, with the following information per sample
        - average coverage
        - PCT coverage >20x
        - Uniformity of coverage (PCT > 0.2*mean) over genome
    """
    coverages = []
    logfiles = glob_files(f"{args.dir}/{sample}/**/{sample}.dragen.bed_coverage_metrics.csv")
    logging.debug(f"List of logfiles to parse: {logfiles}")

    for file in logfiles:
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


def count_cnv(sample):
    """
    Count number of CNVs, by counting line in the VCF file
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame
    """
    cnvs = []
    cnv_dirs = glob_files(f"{sample}.dragen.cnv.vcf.gz")
    logging.debug(f"Files in 'cnv_dirs': {cnv_dirs}")
    count = 0
    for vcf in cnv_dirs:
        path_parts = os.path.split(vcf)
        version    = os.path.basename(path_parts[0])
        with gzip.open(vcf, 'rb') as gzfh:
            gzlines = gzfh.readlines()
        for line in gzlines:
            if line.startswith(b'chr'):
                count += 1
        cnvs.append([sample, version, count])
    return pd.DataFrame(cnvs, columns=['Sample', 'Log NumOfCNVs', 'NumOfCNVs'])


def get_NumOfReads(sample):
    """
    Get Number of Reads for `sample`
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame. There may be multiple NumOfReads files.
                [[Sample, Log, NumOfReads], [], ...]
    """
    files = glob_files(f"{sample}.txt")
    NumOfReads = []
    for file in files:
        path_parts = os.path.split(file)
        dir_parts  = os.path.split(path_parts[0])
        log_NumOfReads = dir_parts[1]
        with open(file, 'r') as fh:
            NumOfReads.append([sample, log_NumOfReads, fh.readline().strip()])
    return pd.DataFrame(NumOfReads, columns=['Sample', 'Log NumOfReads', 'NumOfReads'])


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
    df_coverages = get_coverage_metrics('')
    df_cnvs      = count_cnv('')
    
    # Process list of samples contained in samples_list file.
    #
    samples = get_samples_list(args.samples)
    total   = len(samples)
    for count, sample in enumerate(samples, start=1):
        logging.info(f"Processing {sample}, {count}/{total}")
        logfiles = download_emg_s3_logs(sample, profile='emedgene', logsdir='emg_logs')
        logging.debug(f"S3 logfiles: {logfiles}")
        logfiles     = glob_files(f"{args.dir}/{sample}_v*_sample.log")
        df_metrics   = pd.concat([df_metrics, get_metrics_from_log(sample, logfiles=logfiles)], ignore_index=True)
        df_coverages = pd.concat([df_coverages, get_coverage_metrics(sample)], ignore_index=True)
        df_cnvs      = pd.concat([df_cnvs, count_cnv(sample)], ignore_index=True)

    df_metrics.drop_duplicates(inplace=True)
    df_coverages.drop_duplicates(inplace=True)
    df_cnvs.drop_duplicates(inplace=True)

    df_metrics.reset_index(inplace=True)
    df_coverages.reset_index(inplace=True)
    df_cnvs.reset_index(inplace=True)

    df = df_metrics.merge(df_coverages, on='Sample', how='outer')
    df = df.merge(df_cnvs, on='Sample', how='outer')
    df.drop(['index','index_x', 'index_y'], axis=1, inplace=True)

    df1 = df[['Log filename', 'Log coverage', 'Log NumOfCNVs', 
        'Sample', 'NumOfReads', 'NumOfSNPs', 'NumOfCNVs',
        'Average coverage', 'Coverage uniformity',
        'CNV Number of amplifications', 'CNV Number of passing amplifications', 
        'CNV Number of deletions', 'CNV Number of passing deletions', 
        'PCT coverage >20x',
        'Uniformity of coverage (PCT > 0.2*mean) over genome', 
        'Uniformity of coverage (PCT > 0.4*mean) over genome', 
        'Percent Autosome Callability', 'Estimated sample contamination']]
    df1 = df1.drop_duplicates()

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

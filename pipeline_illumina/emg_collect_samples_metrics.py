#!/usr/bin/env python3
"""
Collect metrics for PRAGMatIQ samples archived on `narval.calculquebec.ca`.

USAGE: emg_get_samples_metrics.py [-d|--directory] 
       emg_get_samples_metrics.py --help

Parse Emedgene's logs files and folders to get analysis metrics for PRAGMatIQ
samples. Download logs from `aws` using the script `archive_PRAGMatIQ.sh`. Log 
files for each sample are archived at: 

`narval.calculquebec.ca:/lustre06/project/6032434/COMMUN/PRAGMatIQ-EMG`

**N.B.** First connect to narval.calculquebec.ca and load environment before 
running this script.

```bash
salloc --account def-rallard --job-name "InteractiveJob" --cpus-per-task 2 --mem-per-cpu 2000  --time 3:0:0
module load scipy-stack/2023a
```
"""

import os
import time
import subprocess
import pandas as pd
from glob import glob

__version__ = 0.1


def parse_args():
    """
    Parse command-line options
    """
    parser = argparse.ArgumentParser(description="Collect metrics for PRAGMatIQ samples archived on `narval.calculquebec.ca`")
    parser.add_argument('-d', '--directory', dest='dir',
                        default="/lustre06/project/6032434/COMMUN/PRAGMatIQ-EMG/archives", 
                        help="Archives directory of all the samples. Default='/lustre06/project/6032434/COMMUN/PRAGMatIQ-EMG/archives'")
    return(parser.parse_args())


def now():
    return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))


def get_NumOfReads(sample):
    """
    Get Number of Reads for `sample`
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame. There may be multiple NumOfReads files.
                [[Sample, Log, NumOfReads], [], ...]
    """
    files = glob(f"{archdir}/{sample}/NumOfReads/*/{sample}.txt")
    NumOfReads = []
    if len(files) > 1:
        print(f"{now()} WARNING: {sample} has multiple NumOfReads")
    for file in files:
        path_parts = os.path.split(file)
        dir_parts  = os.path.split(path_parts[0])
        log_NumOfReads = dir_parts[1]
        with open(file, 'r') as fh:
            NumOfReads.append([sample, log_NumOfReads, fh.readline().strip()])
    return pd.DataFrame(NumOfReads, columns=['Sample', 'Log NumOfReads', 'NumOfReads'])


def get_metrics_from_log(sample):
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

    logfiles = glob(f"{archdir}/{sample}/{sample}_vlocal_*_sample.log")

    # There may be multiple logfiles. Use the filename's _vlocal_ stamp for id
    #
    if len(logfiles) > 1:
        print(f"{now()} WARNING: More than one log file found for {sample}")
    for log in logfiles:
        logname = os.path.basename(log)
        with open(log, "r") as fh:
            lines = fh.readlines()
        reads = ''
        snps  = ''
        cnv_avg_coverage    = ''
        coverage_uniformity = ''
        callability         = ''
        contamination       = '' 
        for line in lines:
            line_parts = line.rstrip().split()

            # For some reason, there are often duplicate lines for these metrics in the log file, 
            # but not in the original bed_coverage_metrics.csv file.
            # The information is sometimes presented with the bytestring (b'blah...\n') symbols 
            # and need to be reformatted.
            #
            if 'Number of reads' in line:
                reads = line_parts[-1].replace("\\n'", "")
            elif 'Average alignment coverage over genome' in line and 'CNV SUMMARY' in line:
                cnv_avg_coverage = line_parts[-1].replace("\\n'", "")
            elif 'Coverage uniformity' in line:
                coverage_uniformity = line_parts[-1].replace("\\n'", "")
            elif 'SNPs' in line and line_parts[11] == "SNPs":
                snps = line_parts[12].replace("\\n'", "")
            elif 'Percent Autosome Callability' in line:
                callability = line_parts[14].replace("\\n'", "")
            elif 'Estimated sample contamination' in line:
                if line_parts[12] != 'standard':
                    contamination = line_parts[12].replace("\\n'", "")
        metrics.append([sample, logname, reads, snps, cnv_avg_coverage, coverage_uniformity, callability, contamination])
    return pd.DataFrame(metrics, columns = ["Sample", "Log filename", "NumOfReads", "NumOfSNPs", "CNV average coverage", "Coverage uniformity", "Percent Autosome Callability", "Estimated sample contamination"])


def count_cnv(sample):
    """
    Count number of CNVs, by counting line in the VCF file
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame
    """
    cnvs = []
    cnv_dirs = glob(f"{archdir}/{sample}/vcf/dragen/*/{sample}.dragen.cnv.vcf.gz")
    count = 0
    for vcf in cnv_dirs:
        path_parts = os.path.split(vcf)
        version    = os.path.basename(path_parts[0])
        vcf_zcat   = subprocess.run(['zcat', vcf], text=True, capture_output=True)
        for line in vcf_zcat.stdout.splitlines():
            if line.startswith('chr'):
                count += 1
        cnvs.append([sample, version, count])
    return pd.DataFrame(cnvs, columns=['Sample', 'Log NumOfCNVs', 'NumOfCNVs'])


def get_coverage_metrics(sample):
    """
    Get coverage metrics for `sample`.
    - `sample`: identifier for sample, ex: "GM231297"
    - Returns : A DataFrame, with the following information per sample
        - average coverage
        - PCT coverage >20x
        - 
    """
    coverages = []
    files = glob(f"{archdir}/{sample}/vcf/dragen/*/{sample}.dragen.bed_coverage_metrics.csv")
    for file in files:
        path_parts   = os.path.split(file)
        version      = os.path.basename(path_parts[0])
        avg_coverage = ''
        coverage_20x = ''
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
                    uniformity_coverage = cols[3]
        coverages.append([sample, version, avg_coverage, coverage_20x, uniformity_coverage])
    return pd.DataFrame(coverages, columns=['Sample', 'Log coverage', 'Average coverage', 'PCT coverage >20x', 'Uniformity of coverage (PCT > 0.2*mean) over genome'])


def main(args):
    """
    From a list of samples, retrieve the following metrics:
    - Institution
    - CQGC_ID
    - Sample_ID
    - #Reads
    - Average Coverage
    - CoverageUniformity
    - %Bases With Coverage > 20X
    - #SNV
    - #Indels
    - #CNVs
    """
    archdir = args.dir
    workdir = os.path.dirname(archdir)
    os.chdir(workdir)


def _test(arg, opt="."):
    print(f"Required command-line argument is: {arg}")
    return(os.stat(opt))


if __name__ == '__main__':
    args = parse_args()
    main(args)
    #_test(args)

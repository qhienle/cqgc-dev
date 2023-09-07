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
                        default="/lustre06/project/6032434/COMMUN/PRAGMatIQ-EMG", 
                        help="Archives directory. Default='/lustre06/project/6032434/COMMUN/PRAGMatIQ-EMG'")
    return(parser.parse_args())


def now():
    return(time.strftime('[%Y-%m-%d@%H:%M:%S]'))


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
    workdir = args.dir
    os.chdir(workdir)


def _test(arg, opt="."):
    print(f"Required command-line argument is: {arg}")
    return(os.stat(opt))


if __name__ == '__main__':
    args = parse_args()
    main(args)
    #_test(args)

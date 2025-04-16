# cqgc-dev

Development of personal utilities for my day job at CQGC. The folder structure of this repo is based on the one for "CQGC-utils.

Workflow: develop herein. When code is ready, copy-pasted individual files to "CQGC-utils", which will be synced with spxp-app02's folder "/staging2/soft/CQGC-utils". 

"./pipeline_illumina/" => "CQGC-utils/Analysis.pipeline_illumina"

## Pre-requisites

```bash
conda deactivate
conda create --name "dev" python=3.13.2
conda activate dev

conda install requests
conda install jupyterlab pandas openpyxl seaborn
conda install -c conda-forge plotly
conda install -c plotly python-kaleido
```
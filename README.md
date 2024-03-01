---
title: "Code for Mathios et al."
author: "AA"
date: "01 March, 2024"
output:
  html_document:
    keep_md: yes
---

# Organization

There are 4 folders of interest in this workflowr.
(1) analysis - contains code needed to generate each figure of the paper.

(2) code - This contains rpcr and rlucas (two packages used extensively in the analysis), some individual r scripts with necessary functions, and 4 folders pertaining to the model creation.

In this paper, we train and evaluate models on 3 cohorts:
1. the LUCAS cohort (158 non-cancer, 129 cancer)
2. the LUCAS cohort excluding prior cancers
3. the LUCAS cohort excluding prior cancers and limited to 50-80 yo smokers.

model_code has the code needed to train models on each cohort, and to perform an external validation of the model from the first cohort. models_c1, models_c2, and models_c3 contain the trained models. For the first two cohorts, two model architectures were trained, for the third cohort only one was trained.

(3) data - this contains raw data used to train the models and generate the figures.
(4) docs - contains html of the markdown files in analysis, as well as the generated figures.

This repository is available on Github, and may be run as a workflowr object to generate a webpage with all code and figures linked.

# Code for figures and tables

All file paths are relative to the top-level directory of this repository.  Below we install non-standard R packages, including supporting packages provided as part of this repository.  Building the SessionInfo file is a useful check that all required packages are available.


```r
install.packages("code/rlucas", repos=NULL, type="source")
install.packages("code/PlasmaTools.lucas", repos=NULL, type="source")
devtools::install_github("jaredhuling/jcolors")
devtools::install_github("/jcolors")
devtools::install_github("cancer-genomics/PlasmaToolsHiseq.hg19")
BiocManager::install(c("paletteer", "precrec"))
wflow_build("analysis/SessionInfo.Rmd")
```

Next, we reproduce the figures and tables in Mathios et al.


```r
wflow_build("analysis/*.Rmd")
```



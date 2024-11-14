LipidIN
===
Introduction
---

Improving annotation accuracy, coverage, speed and depth of lipid profiles remains a significant challenge in traditional spectral matching-based lipidomics. We introduce LipidIN, an advanced framework designed for comprehensive lipid annotation and reverse lipidomics. LipidIN features 168.6 million lipid fragmentation hierarchical library that encompass all potential chain compositions and carbon-carbon double bond locations. Developed expeditious querying module speeds up to around 70 billion times’ spectral querying in less than 1 second. Furthermore, we leverage three relative retention time rules to develop lipid categories intelligence model for reducing false positive annotations and predicting unannotated lipids with a 5.7% estimated false discovery rate coverage 8923 lipids cross various species. More importantly, LipidIN integrates a Wide-spectrum Modeling Yield network for regenerating lipid fingerprints to further improve coverage and accuracy with a 20% estimated recall boosting. The application of LipidIN in multiple tasks demonstrated reliability and potential for lipid annotation and biomarker discovery. 

Demo Watch
---
LipidIN demo watch

[![Video Thumbnail]()](https://github-production-user-asset-6210df.s3.amazonaws.com/148458997/386160788-fd8f9a46-08b2-44fb-8afd-51befcc65f87.mp4?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVCODYLSA53PQK4ZA%2F20241114%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20241114T104854Z&X-Amz-Expires=300&X-Amz-Signature=4fee7de8e430a63f9ea18dddf9cf55af5a51b25ffba928f450ef49914a702826&X-Amz-SignedHeaders=host)

LipidIN(WMYn) demo watch

[![Video Thumbnail]()](https://github-production-user-asset-6210df.s3.amazonaws.com/148458997/386181277-4b9c2441-cc65-49ac-a7a7-5f7d6e47381d.mp4?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVCODYLSA53PQK4ZA%2F20241114%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20241114T115256Z&X-Amz-Expires=300&X-Amz-Signature=122f58bd74e5392ba351627badf5bf41b53cdde19c1dab82b38816e41456aa77&X-Amz-SignedHeaders=host)

General Updates
---
```
Last updated: November 14, 2024

Fixed: (1) Resolved known issues.
Updated: (1) New added 45 kinds of Lipids.
```


```
Last updated: October 08, 2024

Fixed: (1) Resolved the issue of missing header group in some OxPE annotations.
```

```
Last updated: September 09, 2024

Updated: (1) New added Oxidized phosphatidylcholine (OxPC, O2xPC),
                       Oxidized phosphatidylglycerol (OxPG, O2xPG),
                       Oxidized phosphatidylinositol (OxPI, O2xPI),
                       Oxidized phosphatidylserine (OxPS, O2xPS) spectra (https://zenodo.org/records/13735639).
```

```
Last updated: September 03, 2024

Fixed: (1) Update the LCI module to resolve the issue of file accessibility in a multitasking environment;
       (2) Address and rectify errors related to the usage of pos_ALL.rda
           ([https://zenodo.org/records/13645234](https://zenodo.org/records/13645234)).
```

```
Last updated: September 01, 2024

Updated: pos_ALL.rda ([10.5281/zenodo.13624125](https://zenodo.org/records/13624125))
Fixed: (1) New added Osidizid Trilisered (OxTG) spectra.
```
```
Last updated: August 31, 2024

Updated: p1.R; example.R
Fixed: (1) Resolved the prolonged processing time issue in the converting *.mzML format to *.rda format.
```

```
Last updated: August 30, 2024

Added: EQ_support.cpp file and new demos.
Updated: LCI.R; p2.R; p2CH3COO.R; p2COOH.R; example.R
Fixed: (1) Resolved the prolonged processing time issue in the EQ module;
       (2) Fixed the accidental file deletion problem in the LCI module when running multiple files simultaneously.
```


Features
---
        
`Mass Spectrometry Peak Processing Module:`Utilizes existing R packages to process mass spectrometry data in mzML format.
        
`Expeditious querying (EQ) Module:`Performs secondary matching with theoretical or real mass spectrometry libraries and normalizes the matching results.

`The Lipid Categories Intelligence (LCI) Module:`Based on the relative position of primary information, it conducts heuristic searches using secondary matching scores as prior information to re-evaluate high-score matches.

`Reverse Lipid Fingerprint Spectrogram Module:`WMYn predict  lipid fingerprint spectrogram using the model  we designed inspired by KAN and Muit-head attention.

System Architecture
---
The system  main development languages being R，python and C++. While R ,python handles backend processes and C++ accelerates the program. The software employs greedy secondary matching algorithms, heuristic search algorithms, and prior information-based spectrum enhancement algorithms for mass spectrometry analysis, outputting the identification results in CSV format.Using the model to generate reverse lipid fingerprint spectrograms, independent of sample matrices, instruments.

Modules Description and Usage Instructions (LipidIN 4-level hierarchical library)
---
We provide the demo which named "demo neg CH3COO" under LipidIN/LipidIN 4-level hierarchical library, due to github file size limitations, you can find the demos as well as the hierarchical library files "neg_ALL.rda" and "pos_ALL.rda" in the webpage Zenodo (https://zenodo.org/records/13350719).

This task involves searching a 4-level hierarchical library, which is efficient in terms of querying. However, the data format conversion process for the LCI module takes approximately 2 minutes.

To start, locate the file "example.R" within the "LipidIN/LipidIN 4-level hierarchical library" directory. Open this file and modify the relevant parameters as needed. Once you've made the necessary changes to the "example.R" file, select the entire code within the file. Then, click the "Run" button to execute the code.
```
##### Parameter Input #####
FN <- 'I:/LipidIN 4-level hierarchical library/demo neg CH3COO'
pt <- 'I:/LipidIN 4-level hierarchical library'
MS2_filter <- 0.10             
ppm1 <- 5                    
ppm2 <- 10                      
ESI <- 'n2' 
# FN: Address of the *.mzML file to be tested.
# pt: Support code (EQ.cpp, LCI.R, etc.) address.
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter*max intensity will be deleted
```

Modules Description and Usage Instructions (LipidIN MS-DIAL published library)
---
We provide the demo which named "demo of LipidIN" under LipidIN/LipidIN MS-DIAL published library, due to github file size limitations, you can find the demos as well as the file format converted MS-DIAL published library files "MS1_MS2_library.rda" in the webpage Zenodo (https://zenodo.org/records/13350719).

```
##### data preprocessing Using 'RaMS' package #####
source(paste(getwd(),'/preprocessing_RaMS.r',sep=''))
env <- new.env()
preprocessing_RaMS(filename,ESI,MS2_filter)
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter*max intensity will be deleted


##### without multithread data preprocessing Using 'RaMS' package #####
source(paste(getwd(),'/preprocessing_RaMS_nomultithread.r',sep=''))
env <- new.env()
preprocessing_RaMS_nomultithread(filename,ESI,MS2_filter)
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter*max intensity will be deleted
```

`Expeditious querying module (EQ) Module:`Matches mass spectrometry peaks with standard libraries using primary and secondary information.

```
##### EQ module annotation #####
load(paste(getwd(),'/MS1_MS2_library.rda',sep=''))
source(paste(getwd(),'/EQ.r',sep=''))
EQ(filename,ppm1,ppm2,ESI)
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.
# ppm1: MS1 m/z tolerance at parts per million (ppm)
# ppm2: MS2 m/z tolerance at parts per million (ppm)
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
```

`The Lipid Categories Intelligence (LCI) Module:`Reassesses high-confidence matches based on primary information relationships.

```
##### LCI FDR removel #####
source(paste(getwd(),'/LCI.r',sep=''))
env <- new.env()
LCI(filename)
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.


##### without multithread LCI FDR removel #####
source(paste(getwd(),'/LCI_nomultithread.r',sep=''))
env <- new.env()
LCI_nomultithread(filename)
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.
```

`Reverse Lipid Fingerprint Spectrogram Module:`WMYn generates the reverse lipid fingerprint spectrograms.

We provide pre-trained weights. If you need to use them, please download the pre-trained weights and use `predict.py`.

```
if __name__ == "__main__":
    data_path =  " " 
    project_folder = Path(" ")
    mode = " "  # "neg" or "pos"
# data_path is your input data path.
# project_folder is pre-trained weights path.
```

If weights are not available, and for your convenience, we recommend using `train.py`.

```
if __name__ == "__main__":
    data_folder = " "
    output_folder = ' '
    batch_process(data_folder, output_folder)

# data_folder is your input data path.
# output_folder is your output data path.
# The input data include a matrix and a vector for training at least.For example AAA.csv(matrix) and AAA_GT.csv(vector).
# If you need batch processing, please name the corresponding files in the following format, for example, `aaa.csv` with `aaa_GT.csv` ; `bbb.csv` with `bbb_GT.csv`.
# We provide the demo directory, meanwhile, You can use more data. 
```
Requirements
---
All benchmark tests were performed on a personal computer with 13th Gen Intel® Core™ i7-13700F × 16- Core Processor, 64 GB memory, and installed with Windows11 operation system , R-4.2.3 and Python v.3.9 including packages XCMS (v 4.2.2), RaMS (v 1.4.0), parallel (v 3.6.2), doParallel (version 1.0.17), Rcpp (1.0.11), tidyverse (v 1.3.0), WGCNA (v 1.7.0-3), statTarget (v 1.34.0), lightgbm (v 4.5.0), pandas (v 2.0.3.), numpy( v 1.23.5), torch (1.13.1+cu116).
MS entropy and Flash entropy download from Github at https://github.com/YuanyueLi/SpectralEntropy, and https://github.com/YuanyueLi/FlashEntropySearch.
LipidMatch was downloaded from https://github.com/GarrettLab-UF/LipidMatch. 
In all testing, we used MS-DIAL version v4.9.221218 and LipidSearch V4.2. All equations for LipidIN for MS/MS notation and fingerprint regenerating are given in the Methods. Code for clinical cohort analysis  can be access at https://github.com/LinShuhaiLAB/LipidIN/clinical cohort analysis.

- Operating System: Windows
- Compilation Environment
```

Packages need in R code.
> library(this.path)
> library(RaMS)
> library(parallel)
> library(doParallel)
> library(Rcpp)
> library(tidyverse)

Requirements need in python code."
  import os
  import glob
  import torch
  import torch.nn as nn
  import torch.optim as optim
  import pandas as pd
  import random
  import numpy as np
```
- Hardware Requirements: (LipidIN 4-level hierarchical library require RAM >= 16GB, LipidIN MS-DIAL published library has no special requires.)
- Estimated Time for Installing Required Packages
```
  The R code for installing packages takes approximately 10 minutes, depending on individual network conditions.
  The python code for installing requirements takes approximately 10 minutes, depending on individual network conditions.
```
- Estimated Time for Each Module
```
  In the R code：
      Batch realization of converting mzML to rda cost time approximately 10 s
      EQ module cost time approximately 12s (<0.1s for querying and 11 s for data preprocssing of LCI)
      LCI module cost time approximately 1.62 mins
```

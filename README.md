LipidIN
===
Introduction
---
Improving annotation accuracy and coverage, along with obtaining detailed structural information, remains a significant challenge in traditional spectral matching-based lipidomics. We introduce LipidIN, an advanced tool designed for comprehensive lipidome annotation and the generation of personalized reverse lipid fingerprint spectrograms. LipidIN features a hierarchical library with 165.6 million lipids and five-level spectral fragmentation trees that encompass all potential chain compositions and double-bond locations. By incorporating the in-house developed expeditious querying module (EQ), the detailed annotation of lipid structural features, such as head groups, chain lengths, and the degree of unsaturation for each fatty acyl or alkyl/1-alkenyl composition, can be completed at 69.65 billion spectra searches in less than 1 second. Furthermore, we introduced three rules for relative retention times: equivalent carbon number (ECN), interclass unsaturation parallelism (IUP), and equivalent separated carbon number (ESCN), to develop the Lipid Categories Intelligence (LCI) model for eliminating false positive annotations and predicting unidentified lipids. LipidIN also integrates a Wide-spectrum Modeling Yield network model (WYMn), generating reverse lipid fingerprint spectrograms, independent of sample matrices, instruments, and analytical methods. LipidIN  providing a robust and efficient resource for lipidomics research. Its reliability has been validated through 1,833 clinical cohorts, underscoring its effectiveness in clinical lipidomics research.

Features
---
        
`Mass Spectrometry Peak Processing Module:`Utilizes existing R packages to process mass spectrometry data in mzML format.
        
`Expeditious querying module (EQ) Module:`Performs secondary matching with theoretical or real mass spectrometry libraries and normalizes the matching results.

`The Lipid Categories Intelligence (LCI) Module:`Based on the relative position of primary information, it conducts heuristic searches using secondary matching scores as prior information to re-evaluate high-score matches.

`Reverse Lipid Fingerprint Spectrogram Module:`WMYn predict  lipid fingerprint spectrogram using the model  we designed inspired by KAN and Muit-head attention.

System Architecture
---
The system  main development languages being R，python and C++. While R ,python handles backend processes and C++ accelerates the program. The software employs greedy secondary matching algorithms, heuristic search algorithms, and prior information-based spectrum enhancement algorithms for mass spectrometry analysis, outputting the identification results in CSV format.Using the model to generate reverse lipid fingerprint spectrograms, independent of sample matrices, instruments.

Modules Description and Usage Instructions
---
![Fig 1-1](https://github.com/user-attachments/assets/7c2470c4-7609-40fe-9ce5-c7167dca8aa9)

We provide the demo which named demo of LipidIN (without .mzML) under LipidIN/MS-DIAL published library 

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

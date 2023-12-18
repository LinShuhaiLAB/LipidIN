# LipidIN (Lipidomics Integration)


**Table of Contents**

[TOCM]

[TOC]



Introduction
=============
This software is an automated lipid profiling tool that aims to help medical and biological laboratories and clinics more rapidly and accurately diagnose a variety of bio-polar metabolites. The software is compatible with various instruments currently available on the market. The software utilizes advanced models, algorithms and artificial intelligence techniques to automatically analyze the composition of various lipids in samples. By using this software, users do not need to establish standard spectra libraries for identification and can save a great deal of time and effort while improving the accuracy and reliability of qualitative analysis. Compared to traditional manual qualitative methods, this new software-based approach significantly enhances workflow efficiency and diagnostic capabilities. Preliminary evaluations indicate that the software is able to reliably identify major lipid subclasses and species from complex biological extracts. Overall, the automated lipid profiling tool shows promise for streamlining lipidomic workflows in clinical diagnostics and biomarker discovery. Further validation and large-scale prospective studies are warranted to fully characterize the performance and utility of this new informatics platform.

Overall Program Workflow
=============
The workflow of the software program can be divided into four main steps:

Users submit files and parameters for preprocessing in the peak processing module to initiate lipid identification requests. The peak processing module receives the request and uses XCMS and Msnbase for data matching according to the submitted files and parameters.

Users upload the theoretical library and parameters for unsupervised secondary matching in the secondary matching module without a priori information. This module performs greedy unsupervised secondary matching on the preprocessed files and classifies the mass spectral peaks according to the matching scores.

Users input parameters in the heuristic search module. This module performs heuristic searches on the secondary matching results according to the user-input parameters and returns the loss function scores for each peak result, finally outputting the possible results for each mass spectral peak and the reference scores for each result.

Users input parameters in the secondary matching module with a priori information. This module predicts the spectra for mass spectral peaks with poor quality according to the user-input parameters and finally outputs the possible results for each mass spectral peak and the reference scores for each result.

Guidelines for Use
=============
To operate the automated lipid profiling software:

Launch the graphical user interface by clicking the "Run App" button in the top right of RStudio.

The Shiny web application interface will launch in a browser window/tab after a short delay.

Use the browse button to select the raw mass spectrometry data file to profile. Supported file formats include mzML, mzXML and fid.

Adjust additional analysis parameters if desired, such as preprocessing and identification settings. Default optimized values are pre-populated.

Click "Submit" to begin the automated profiling workflow, including preprocessing, feature extraction, database matching and lipid identification.

Identification results will be displayed in a sortable table, showing all putative matches along with associated confidence scores.

Individual spectra can be visualized alongside assignments for quality control.

Save or export results files as needed in common formats like CSV or TXT.

Troubleshoot any issues by referring to the documentation accessible via the help menu. Contact technical support for assistance as needed.

By following these simple steps, users can quickly generate lipid identifications from their raw MS data without specialist bioinformatics skills.










### End

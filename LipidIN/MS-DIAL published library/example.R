# Preparation and installation packages ----------------------------------------
packages <- c('this.path','parallel','doParallel','RaMS','Rcpp','tidyverse')
installed_packages <- packages %in% rownames(installed.packages())
if(all(installed_packages)){
  print("All required packages are installed.")
} else{
  print("The following packages are missing:")
  print(packages[!installed_packages])
  
  # Install the missing packages
  install.packages(packages[!installed_packages])
}
library(this.path)
library(RaMS)
library(parallel)
library(doParallel)
library(Rcpp)
library(tidyverse)
##### data preprocessing Using 'RaMS' package #####
source(paste(getwd(),'/preprocessing_RaMS.r',sep=''))
env <- new.env()
# example
a <- Sys.time()
filename <- 'D:/Draw/LipidIN MS-DIAL published library/demo of LipidIN/QCL_NEG_ID_1.mzML'
ESI <- 'n2'
MS2_filter <- 0.10
preprocessing_RaMS(filename,ESI,MS2_filter)
Sys.time()-a
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter*max intensity will be deleted


##### without multithread data preprocessing Using 'RaMS' package #####
source(paste(getwd(),'/preprocessing_RaMS_nomultithread.r',sep=''))
env <- new.env()
# example
a <- Sys.time()
filename <- 'D:/Draw/LipidIN MS-DIAL published library/demo of LipidIN/QCL_NEG_ID_1.mzML'
ESI <- 'n2'
MS2_filter <- 0.10
preprocessing_RaMS_nomultithread(filename,ESI,MS2_filter)
Sys.time()-a
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter*max intensity will be deleted


##### EQ module annotation #####
load(paste(getwd(),'/MS1_MS2_library.rda',sep=''))
# example
a <- Sys.time()
source(paste(getwd(),'/EQ.r',sep=''))
filename <- 'D:/Draw/LipidIN MS-DIAL published library/demo of LipidIN/QCL_NEG_ID_1.mzML'
ppm1 <- 5
ppm2 <- 10
ESI <- 'n2'
EQ(filename,ppm1,ppm2,ESI)
Sys.time()-a
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.
# ppm1: MS1 m/z tolerance at parts per million (ppm)
# ppm2: MS2 m/z tolerance at parts per million (ppm)
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.


##### LCI FDR removel #####
source(paste(getwd(),'/LCI.r',sep=''))
env <- new.env()
# example
a <- Sys.time()
filename <- 'D:/Draw/LipidIN MS-DIAL published library/demo of LipidIN/QCL_NEG_ID_1.mzML'
LCI(filename)
Sys.time()-a
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.


##### without multithread LCI FDR removel #####
source(paste(getwd(),'/LCI_nomultithread.r',sep=''))
env <- new.env()
# example
a <- Sys.time()
filename <- 'D:/Draw/LipidIN MS-DIAL published library/demo of LipidIN/QCL_NEG_ID_1.mzML'
LCI_nomultithread(filename)
Sys.time()-a
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.











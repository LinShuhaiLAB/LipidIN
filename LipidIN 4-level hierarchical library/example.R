FN <- 'D:/bio_inf/LipidIN/LipidIN 4-level hierarchical library v1.0.4/demo pos'
pt <- 'D:/bio_inf/LipidIN/LipidIN 4-level hierarchical library v1.0.4'
MS2_filter <- 0.10             
ppm1 <- 5                    
ppm2 <- 10                      
ESI <- 'p' 
# FN: Address of the *.mzML file to be tested.
# pt: Support code (EQ.cpp, LCI.R, etc.) address.
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter*max intensity will be deleted

# Preparation and installation packages ----------------------------------------
packages <- c('this.path','parallel','doParallel','RaMS','Rcpp','tidyverse','dplyr')
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
library(dplyr)
source(paste(pt,'/p1.r',sep=''))
sourceCpp(paste(pt,'/EQ.cpp',sep=''))
sourceCpp(paste(pt,'/EQ_support.cpp',sep=''))
sourceCpp(paste(pt,'/calculateMolecularMass.cpp',sep=''))
sourceCpp(paste(pt,'/rawmzcluster.cpp',sep=''))
sourceCpp(paste(pt,'/removeRowsWithinError.cpp',sep=''))
sourceCpp(paste(pt,'/sortMatrixByRow.cpp',sep=''))

# Batch realization of converting mzML to rda ------------------------------
setwd(FN)
FN1 <- list.files()
FN1 <- FN1[grep('.mzML',FN1)]
a <- Sys.time()
env <- new.env()
for(ii in FN1){
  print(ii)
  part1_RaMS(ii,ESImode=ESI,MS2_filter)
}
Sys.time()-a
# EQ module --------------------------------------------------------------------
if(ESI=='p'){
  load(paste(pt,'/pos_ALL.rda',sep=''))
  source(paste(pt,'/p2.r',sep=''))
}
if(ESI!='p'){
  load(paste(pt,'/neg_ALL.rda',sep=''))
  if(ESI=='n1'){
    source(paste(pt,'/p2COOH.r',sep=''))
  }
  if(ESI=='n2'){
    source(paste(pt,'/p2CH3COO.r',sep=''))
  }
}
FN1 <- list.files()
FN1 <- FN1[grep('.rda',FN1)]
a <- Sys.time()
compare <- new(Comparator)
compare$Load(NULL,count.data.frame$MZ,ppm1,ppm2,count.list)
for(ii in FN1){
  print(ii)
  part2(ii)
}
Sys.time()-a

# LCI module -------------------------------------------------------------------
setwd(pt)
source('LCI.r')
rm(count.data.frame,count.list)
a <- Sys.time()
env <- new.env()
for(ii in FN1){
  print(ii)
  LCI(paste(FN,'/',ii,sep=''))
}
Sys.time()-a
setwd(FN)
a <- list.files()
b <- a[grep('.csv',a)]
b <- b[-grep('final_output.csv',b)]
file.remove(b)













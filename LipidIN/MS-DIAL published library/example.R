##### data preprocessing Using 'RaMS' package #####
env=new.env()
source(paste(getwd(),'/preprocessing_RaMS.r',sep=''))
preprocessing_RaMS(filename,ESI,MS2_filter)
# filename: Location of .mzML file, for example '.../demo pos/QC_POS1.mzML'.
# ESI: 'p' for positive ionization mode，
#      ‘n1’ for negative ionization mode [M+COOH]-，
#      'n2' for negative ionization mode [M+CH3COO]-.
# MS2_filter: a value of 0-1, MS2 fragments with intensity lower than the MS2_filter will be deleted

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

# LCI去除假阳性 ----------------------------------------------------------------
source(paste(getwd(),'/LCI.r',sep=''))
env <- new.env()
LCI(filename)
# filename: Location of .rda file output by data preprocessing, for example '.../demo pos/QC_POS1.rda'.












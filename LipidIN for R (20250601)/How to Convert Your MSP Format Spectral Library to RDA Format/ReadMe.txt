What you need to do is open the R script file (How to convert your file format.R) and modify the following sections according to your MSP file.
# Modify your parameters here ---------------------------------------------
# install.packages("pbapply")
FN1 <- 'D:/xxxxxxxx'
# The address where the working directory is located.
FN2 <- 'MoNA-export-All_LC-MS-MS_Orbitrap.msp'
# The name of the MSP file for the spectral library.
Meta_name <- 'Name: '
# Metabolite or lipid name prompt.
Formula_label <- 'Formula: '
# Formula prompt.
PrecursorMZ_label <- 'PrecursorMZ: '
# Precursor Ion prompt.
Precursor_type_label <- 'Precursor_type'
# Amalgamation prompt.
Ion_mode_label <- 'Ion_mode: '
# Ionization mode prompt.
Num_Peaks_label <- 'Num Peaks: '
# Prompt for the number of secondary peaks.
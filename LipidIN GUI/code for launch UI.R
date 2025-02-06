packages <- c("this.path", 'magick',"shinyjs", "RaMS", "Rcpp", "tidyverse", "dplyr", "shiny","golem", "shinydashboard", "ggplot2", "ggforce", "magick", "grid",'devtools',"tidyr", "shinyFiles")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
    message(paste(pkg, "has been installed"))
  } else {
    message(paste(pkg, "is already installed"))
  }
}


library(devtools)

install.packages(
  pkgs = 'D:/bio_inf/LipidIN-main/LipidIN GUI/LipidIN_2.0.0.2.tar.gz',
  # Note that here you need to write the full address 
  # of LipidIN_2.0.0.1.tar.gz, the above is only an example, 
  # please modify the address according to the actual situation.
  lib = .libPaths()[length(.libPaths())],
  repos = NULL,
  dependencies = T
)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("xcms", quietly = TRUE))
  BiocManager::install("xcms")

if (!require("CAMERA", quietly = TRUE))
  BiocManager::install("CAMERA")

library(LipidIN)
run_LipidIN()

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
  pkgs = '~/LipidIN-main/LipidIN GUI/LipidIN_2.0.0.1.tar.gz',
  lib = .libPaths()[length(.libPaths())],
  repos = NULL,
  dependencies = T
)


library(LipidIN)
run_LipidIN()

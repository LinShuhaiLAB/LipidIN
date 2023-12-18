# install packages
pkg <- c('propagate',
         'investr',
         'tidyr',
         'tidyverse',
         'dplyr',
         'Rcpp',
         'shiny',
         'shinydashboard',
         'ggplot2',
         'ggrepel',
         'propagate',
         'investr',
         'ggrepel'
         )
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("MSnbase", quietly = TRUE))
  BiocManager::install("MSnbase")
if (!require("xcms", quietly = TRUE))
  BiocManager::install("xcms")
if (!require("bs4Dash", quietly = TRUE))
  BiocManager::install("bs4Dash")
setwd('E:/metaID/20230810newui/function')
library(propagate)
library(investr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(MSnbase)
library(xcms)
library(Rcpp)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggrepel)
library(propagate)
library(investr)
library(ggrepel)
source('sourceFunction.R')
sourceCpp('removeRowsWithinError.cpp')
sourceCpp('all.cpp')
sourceCpp("fastMatch.cpp")
sourceCpp('calculateMolecularMass.cpp')
sourceCpp('sortMatrixByRow.cpp')
source('ui.r')
source('server.r')

# 启动本地web端
shinyApp(ui=ui,server=server)













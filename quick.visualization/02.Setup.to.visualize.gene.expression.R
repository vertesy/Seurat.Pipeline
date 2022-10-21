######################################################################
# SETUP for Visualize single-cell gene expression
######################################################################
# source ('')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)




# Load packages ------------------------
require(tidyverse) # graphics / utils
require(cowplot)
require(Seurat)
require(tictoc)
# require(schex)
require(clipr)


# 2. Load custom function libraries ------------------------
require(Stringendo)
require(ReadWriter)
require(CodeAndRoll2) # source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R')
require(MarkdownHelpers)
require(MarkdownReports)
require(ggExpress);
require(Seurat.utils)

# 3. Source another custom library that is not a package ------------------------
{
  Rocinante.local <- '~/GitHub/Packages/Rocinante/R/Rocinante.R'
  Rocinante.https <- 'https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R'
  if (file.exists(Rocinante.local)) {
    print('Rocinante.local')
    try(source(Rocinante.local), silent = F)
  } else if (RCurl::url.exists(Rocinante.https)) {
    print('Rocinante.https')
    try(source('https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R'), silent = F)
  } else {
    print('ROCINANTE NOT FOUND!!!')
  }
}


# ------------------------------------------------------------------------
print("All packages loaded")




# Setup & Custom functions you need--------------------------------------------------------------------------------
source("https://raw.githubusercontent.com/vertesy/TheCorvinas/master/R/CodeAndRoll.R")
source("https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/Seurat3.Functions.R")


# require(cowplot)
require(tidyverse)
require(Seurat)
require(ggplot2)
require(MarkdownReportsDev)
  # to install MarkdownReportsDev:
  # require(devtools)
  # devtools::install_github(repo = "vertesy/MarkdownReportsDev")



# Metadata ------------------------
# Parameters ------------------------



# Install ------------------------

if (FALSE) {
  source("https://raw.githubusercontent.com/vertesy/TheCorvinas/master/R/CodeAndRoll.R")
  source("https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/Seurat3.Functions.R")

  install.packages("withr"); # If you don't have it
  require("devtools")
  require(withr)
  # devtools::install_github(repo = "vertesy/MarkdownReportsDev")
  with_libpaths(new = "C:/Program Files/R/R-3.6.1/library", install_github(repo = "vertesy/MarkdownReportsDev"))

  install.packages("tidyverse")
  install.packages("ggplot2")
  install.packages('tictoc')
  install.packages('Seurat')
  install.packages('gtools')
  install.packages('clipr')
  install.packages('doMC')
  install.packages('readr')
  # install.packages("gdata")
  # install.packages(doMC)
  # install.packages("doMC", repos = "http://R-Forge.R-project.org", lib = "C:/Program Files/R/R-3.6.1/library" )

}


######################################################################
# Load.packages.local.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/Load.packages.local.R')
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/Load.packages.local.R')


# 1. Load packages ------------------------
require(tidyverse) # graphics / utils
require(cowplot)
require(colorout)

require(Seurat)

require(future) # parallelization
# require(doMC)

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
require(Seurat.utils) # try(source("~/GitHub/Packages/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R"));
# try(source("~/GitHub/Packages/Seurat.multicore/Seurat3.Multicore.Load.R"));
require(UVI.tools)

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






# source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R'),silent = F)
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/00.Load.Seurat.Utils.WEB.R'), silent =   T)
# try(source('https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R?token=ABG5SV5KJJH7L7TJS66245K7Y65QQ'), silent =   F)
## try(source('https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/00.Load.Seurat3.Multicore.WEB.R'), silent=T)


# source('https://raw.githubusercontent.com/vertesy/DatabaseLinke.R/master/DatabaseLinke.R')

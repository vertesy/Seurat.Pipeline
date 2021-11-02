######################################################################
# Load.packages.CBE.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Load.packages.CBE.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.CBE.R")

# 1. Load packages ------------------------
require(tidyverse) # graphics / utils
require(cowplot)
require(MarkdownReports)
require(colorout)

require(Seurat)

require(future) # parallelization
require(doMC)

require(tictoc)
# require(schex)

# 2. Load custom function libraries ------------------------
try(require(CodeAndRoll2)
require(Stringendo)
require(ReadWriter)
# source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R'),silent = F)
# require(CodeAndRoll2)
require(Stringendo)
require(ReadWriter)
# source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R')
try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/00.Load.Seurat.Utils.WEB.R'), silent =   T)
# try(source("~/GitHub/Packages/Seurat.utils/00.Load.Seurat.Utils.WEB.R"));
try(source('https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R'), silent =   F)
# require(ggExpress);

# try(source("~/GitHub/Packages/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R"));
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/00.Load.Seurat3.Multicore.WEB.R'), silent=T)
# # try(source("~/GitHub/Packages/Seurat.multicore/Seurat3.Multicore.Load.R"));

# 3. Load CBE specific function libraries ------------------------
if (IfExistsAndTrue("onCBE")) {
  try(dyn.load("/software/2020/software/nodejs/10.15.1-foss-2018b/lib/libnode.so.64"), silent = F)
  try(library(V8,lib.loc = '/software/2020/software/v8/3.3.1-foss-2018b-r-4.0.2/'), silent = F)
  require("htmltools")

  # options for CBE ------------------------
  options(bitmapType = 'cairo')

} else {
  print("!!! onCBE is not defined as TRUE, thus CBE specific settings are not loaded.")
}



print("Packages loaded")


######################################################################
# Load.packages.CBE.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/Load.packages.CBE.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/Load.packages.CBE.R")

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

source('https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R')

if (!require(Stringendo)) source('https://raw.githubusercontent.com/vertesy/Stringendo/main/R/Stringendo.R');
if (!require(ReadWriter)) source('https://raw.githubusercontent.com/vertesy/ReadWriter/main/R/ReadWriter.R');
if (!require(CodeAndRoll2)) source('https://raw.githubusercontent.com/vertesy/CodeAndRoll2/main/R/CodeAndRoll2.R');
if (!require(MarkdownHelpers)) source('https://raw.githubusercontent.com/vertesy/MarkdownHelpers/main/R/MarkdownHelpers.R');
if (!require(MarkdownReports)) source('https://raw.githubusercontent.com/vertesy/MarkdownReports/master/R/MarkdownReports.R');
if (!require(ggExpress)) source('https://raw.githubusercontent.com/vertesy/ggExpress/master/R/ggExpress.R');
if (!require(Seurat.utils)) source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/R/Seurat.Utils.R');
if (!require(UVI.tools)) source('https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.R?token=GHSAT0AAAAAACHD3DED5LL7SCKYQ2LEEXGIZHVWHCQ');
if (!require(UVI.tools)) source('https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.Bulk.R?token=GHSAT0AAAAAACHD3DEDZA3IEHUTC2ICHB7EZHVWHEA');
if (!require(Connectome.tools)) source('https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.R?token=GHSAT0AAAAAACHD3DEDG2MJ54KUDV6BIPESZHVWDNQ');
if (!require(Connectome.tools)) source('https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.AAV.R?token=GHSAT0AAAAAACHD3DEDDISG3JIQWPSNADBKZHVWDWQ');


# 3. Load CBE specific function libraries ------------------------
if (IfExistsAndTrue("onCBE")) {
  "not working, pr"
  try(dyn.load("/software/2020/software/nodejs/10.15.1-foss-2018b/lib/libnode.so.64"), silent = F)
  try(library(V8,lib.loc = '/software/2020/software/v8/3.3.1-foss-2018b-r-4.0.2/'), silent = F)
  require("htmltools")

  # options for CBE ------------------------
  options(bitmapType = 'cairo')

} else {
  print("!!! onCBE is not defined as TRUE, thus CBE specific settings are not loaded.")
}



print("Packages loaded")


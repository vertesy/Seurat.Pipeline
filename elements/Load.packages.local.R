######################################################################
# Load.packages.local.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Load.packages.local.R')
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.local.R')

# Load packages ------------------------
require(tidyverse) # graphics / utils
require(cowplot)
require(MarkdownReportsDev)
require(colorout)

require(Seurat)

require(future) # parallelization
require(doMC)

require(tictoc)
require(schex)

# 2. Load function libraries ------------------------
try(source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'))
try(source("~/GitHub/Packages/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R"));
try(source("~/GitHub/Packages/ggExpressDev/ggExpress.functions.R"));
# try(source("~/GitHub/Packages/Seurat.multicore/Seurat3.Multicore.Load.R"));

# try(source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R'),silent = F)
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/00.Load.Seurat.Utils.WEB.R'), silent =   T)
# try(source('https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R?token=ABG5SV5KJJH7L7TJS66245K7Y65QQ'), silent =   F)
## try(source('https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/00.Load.Seurat3.Multicore.WEB.R'), silent=T)

print("Packages loaded")


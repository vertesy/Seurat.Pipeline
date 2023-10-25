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

load_or_source <- function(path, web) {
  if (dir.exists(path)) devtools::load_all(path) else {
    print('Local package not found'); print(path); print('Sourcing from web')
    source(web)
  }
}

# 2. Load custom function libraries ------------------------
source('~/GitHub/Packages/Rocinante/R/Rocinante.R')
# source('https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R')

load_or_source(path = '~/GitHub/Packages/Stringendo', web =  'https://raw.githubusercontent.com/vertesy/Stringendo/main/R/Stringendo.R'); warnings()
load_or_source(path = '~/GitHub/Packages/ReadWriter', web =  'https://raw.githubusercontent.com/vertesy/ReadWriter/main/R/ReadWriter.R'); warnings()
load_or_source(path = '~/GitHub/Packages/CodeAndRoll2', web =  'https://raw.githubusercontent.com/vertesy/CodeAndRoll2/main/R/CodeAndRoll2.R'); warnings()
load_or_source(path = '~/GitHub/Packages/MarkdownHelpers', web =  'https://raw.githubusercontent.com/vertesy/MarkdownHelpers/main/R/MarkdownHelpers.R'); warnings()
load_or_source(path = '~/GitHub/Packages/MarkdownReports', web =  'https://raw.githubusercontent.com/vertesy/MarkdownReports/master/R/MarkdownReports.R'); warnings()
load_or_source(path = '~/GitHub/Packages/ggExpress', web =  'https://raw.githubusercontent.com/vertesy/ggExpress/master/R/ggExpress.R'); warnings()
load_or_source(path = '~/GitHub/Packages/Seurat.utils', web =  'https://raw.githubusercontent.com/vertesy/Seurat.utils/master/R/Seurat.Utils.R'); warnings()

"Below ones wont work"
load_or_source(path = '~/GitHub/Packages/UVI.tools', web =  'https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.R?token=    ');
# load_or_source(path = '~/GitHub/Packages/UVI.tools', web =  'https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.Bulk.R?token=    ');
load_or_source(path = '~/GitHub/Packages/Connectome.tools', web =  'https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.R?token=    ');
# load_or_source(path = '~/GitHub/Packages/Connectome.tools', web =  'https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.AAV.R?token=    ');

# Old way ------------------------

# if (!require(Stringendo)) source('https://raw.githubusercontent.com/vertesy/Stringendo/main/R/Stringendo.R');
# if (!require(ReadWriter)) source('https://raw.githubusercontent.com/vertesy/ReadWriter/main/R/ReadWriter.R');
# if (!require(CodeAndRoll2)) source('https://raw.githubusercontent.com/vertesy/CodeAndRoll2/main/R/CodeAndRoll2.R');
# if (!require(MarkdownHelpers)) source('https://raw.githubusercontent.com/vertesy/MarkdownHelpers/main/R/MarkdownHelpers.R');
# if (!require(MarkdownReports)) source('https://raw.githubusercontent.com/vertesy/MarkdownReports/master/R/MarkdownReports.R');
# if (!require(ggExpress)) source('https://raw.githubusercontent.com/vertesy/ggExpress/master/R/ggExpress.R');
# if (!require(Seurat.utils)) source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/R/Seurat.Utils.R');
# if (!require(UVI.tools)) source('https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.R?token=    ');
# if (!require(UVI.tools)) source('https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.Bulk.R?token=    ');
# if (!require(Connectome.tools)) source('https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.R?token=    ');
# if (!require(Connectome.tools)) source('https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.AAV.R?token=    ');


# 3. Load CBE specific function libraries ------------------------
if (IfExistsAndTrue("onCBE")) {
  # "not working, pr"
  # try(dyn.load("/software/2020/software/nodejs/10.15.1-foss-2018b/lib/libnode.so.64"), silent = F)
  # try(library(V8,lib.loc = '/software/2020/software/v8/3.3.1-foss-2018b-r-4.0.2/'), silent = F)
  # require("htmltools")

  # options for CBE ------------------------
  options(bitmapType = 'cairo')
  oo = list.files
  b.raster = TRUE
} else {
  print("!!! onCBE is not defined as TRUE, thus CBE specific settings are not loaded.")
}


print("Custom packages loaded!")




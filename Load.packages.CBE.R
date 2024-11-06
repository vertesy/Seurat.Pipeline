######################################################################
# Load.packages.CBE.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/Load.packages.CBE.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/Load.packages.CBE.R")
# source('~/.pack.R')

onCBE = TRUE; # stopifnot(exists('onCBE'))

# 1. Load packages ------------------------
require(tidyverse) # graphics / utils
require(cowplot)
require(colorout)

require(Seurat)

require(future) # parallelization
require(doMC)
require(checkmate)

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
# load_or_source(path = '~/GitHub/Packages/gruffi', web =  'https://raw.githubusercontent.com/jn-goe/gruffi/main/R/gruffi.R'); warnings()
# devtools::install_github(repo = "jn-goe/gruffi", upgrade = F)

load_or_source(path = '~/GitHub/Packages/isoENV', web =  'https://raw.githubusercontent.com/vertesy/isoENV/main/R/isoENV.R');

"Below ones wont work online"
load_or_source(path = '~/GitHub/Packages/NestedMultiplexer', web =  'https://raw.githubusercontent.com/vertesy/NestedMultiplexer/main/R/NestedMultiplexer?token=    ');
load_or_source(path = '~/GitHub/Packages/UVI.tools', web =  'https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.R?token=    ');
# load_or_source(path = '~/GitHub/Packages/UVI.tools', web =  'https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.Bulk.R?token=    ');
load_or_source(path = '~/GitHub/Packages/Connectome.tools', web =  'https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.R?token=    ');
# load_or_source(path = '~/GitHub/Packages/Connectome.tools', web =  'https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.AAV.R?token=    ');

# try(load_or_source(path = '~/GitHub/Packages/PackageTools', web =  'https://raw.githubusercontent.com/vertesy/PackageTools/main/R/PackageTools.R'), silent = T)

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



print("Custom packages loaded!")




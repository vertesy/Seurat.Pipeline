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
require(devtools)


# 2. Load custom function libraries ------------------------

load_all(path = '~/GitHub/Packages/Stringendo');
load_all(path = '~/GitHub/Packages/ReadWriter');
load_all(path = '~/GitHub/Packages/CodeAndRoll2');
load_all(path = '~/GitHub/Packages/MarkdownHelpers');
load_all(path = '~/GitHub/Packages/MarkdownReports');
load_all(path = '~/GitHub/Packages/ggExpress');
load_all(path = '~/GitHub/Packages/Seurat.utils');

load_all(path = '~/GitHub/Packages/UVI.tools');
load_all(path = '~/GitHub/Packages/Connectome.tools');


#
# require(Stringendo) # source('https://raw.githubusercontent.com/vertesy/Stringendo/main/R/Stringendo.R')
# require(ReadWriter) # source('https://raw.githubusercontent.com/vertesy/ReadWriter/main/R/ReadWriter.R')
# require(CodeAndRoll2) # source('https://raw.githubusercontent.com/vertesy/CodeAndRoll2/main/R/CodeAndRoll2.R')
# require(MarkdownHelpers) # source('https://raw.githubusercontent.com/vertesy/MarkdownHelpers/main/R/MarkdownHelpers.R')
# require(MarkdownReports) # source('https://raw.githubusercontent.com/vertesy/MarkdownReports/master/R/MarkdownReports.R')
# require(ggExpress) # source('https://raw.githubusercontent.com/vertesy/ggExpress/master/R/ggExpress.R')
# require(Seurat.utils) # source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/R/Seurat.Utils.R')
# require(UVI.tools) # source('https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.R?token=GHSAT0AAAAAACHD3DED5LL7SCKYQ2LEEXGIZHVWHCQ')
# # source('https://raw.githubusercontent.com/vertesy/UVI.tools/main/R/UVI.tools.Bulk.R?token=GHSAT0AAAAAACHD3DEDZA3IEHUTC2ICHB7EZHVWHEA')
# require(Connectome.tools) # source('https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.R?token=GHSAT0AAAAAACHD3DEDG2MJ54KUDV6BIPESZHVWDNQ')
# # source('https://raw.githubusercontent.com/vertesy/Connectome.tools/main/R/Connectome.tools.AAV.R?token=GHSAT0AAAAAACHD3DEDDISG3JIQWPSNADBKZHVWDWQ')

# source('https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R')

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

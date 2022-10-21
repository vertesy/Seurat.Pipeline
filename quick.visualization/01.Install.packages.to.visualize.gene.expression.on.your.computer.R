######################################################################
# Install packages to Visualize single-cell gene expression
# Fri Oct 21 10:33:06 2022
######################################################################
# source ('')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)




# Install ------------------------


install.packages("tidyverse")
install.packages('tictoc')
install.packages('Seurat')
install.packages('gtools')
install.packages('clipr')
install.packages('readr')
install.packages("devtools")
install.packages('BiocManager')
install.packages('cowplot')


# Install from Bioconductor packages ------------------------
BiocManager::install("sparseMatrixStats")

# Install custom packages ------------------------
require("devtools")

# Install dependencies
devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)

# Recommended
devtools::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)

# Install Seurat.utils
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)


"If you see errors that you miss any additional package, you can install it just like the ones above."

# Source another custom library that is not a package ------------------------
print('sourcing Rocinante from the web')
try(source('https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R'), silent = F)



######################################################################
# plot.3D.umaps.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/plot.3D.umaps.R")

# Functions ------------------------
library(plotly)
source('~/GitHub/Packages/Seurat.utils/plotting.dim.reduction.3D.R')
create_set_Original_OutDir()



# Setup ------------------------

combined.obj <- RecallReduction(obj = combined.obj, reduction = "umap", dim = 3)
(categ3Dplots <- unique(c(names(tags), "Phase", GetClusteringRuns(),  GetNamedClusteringRuns() )))


# Plot Categories ------------------------
GetNamedClusteringRuns()
Plot3D.ListOfCategories(obj = combined.obj, cex = 1.75, ListOfCategories = categ3Dplots, annotate.by =  GetNamedClusteringRuns())


# Plot Genes ------------------------
genes.for.3D.plots <- union.ls(list(ClassicMarkers.plus, combined.obj@misc[[ppp('top.markers.res', p$def_res)]]))
combined.obj <- AddGOGeneList.manual(genes = genes.for.3D.plots, GO = "genes.for.3D.plots")
Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = genes.for.3D.plots, annotate.by =  GetNamedClusteringRuns()[1], cex = 2)


# End ------------------------
combined.obj <- RecallReduction(obj = combined.obj, reduction = "umap", dim = 2)
OutDir = getwd()

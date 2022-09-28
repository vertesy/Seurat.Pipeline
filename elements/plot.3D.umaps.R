######################################################################
# plot.3D.umaps.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/plot.3D.umaps.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Plot.3D.umaps.R")

# Functions ------------------------
library(plotly)
# source('~/GitHub/Packages/Seurat.utils/Functions/Plotting.dim.reduction.3D.R')
create_set_Original_OutDir()



# Setup ------------------------

combined.obj <- RecallReduction(obj = combined.obj, reduction = "umap", dim = 3)




# Plot Classic Marker Genes ------------------------
if (exists('genes.ls')) {
  ClassicMarkers.for.3D.plots <- union.ls(list(genes.ls$'ClassicMarkers'))
  combined.obj <- gruffi::AddGOGeneList.manual(genes = ClassicMarkers.for.3D.plots, GO = "ClassicMarkers.for.3D.plots")
  Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = ClassicMarkers.for.3D.plots, annotate.by =  GetNamedClusteringRuns()[1], cex = 2)
}

# Plot DE Genes ------------------------
top.markers.exist <- any(grepl('top.markers.res.', names(combined.obj@misc)))
if (top.markers.exist) {
  DE.genes.for.3D.plots <- combined.obj@misc[[ppp('top.markers.res', p$"def_res")]]
  combined.obj <- gruffi::AddGOGeneList.manual(genes = DE.genes.for.3D.plots, GO = "DE.genes.for.3D.plots")
  Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = DE.genes.for.3D.plots, annotate.by =  GetNamedClusteringRuns()[1], cex = 2)
}


# Plot Categories ------------------------
{
  if(!exists("meta.tags")) list.fromNames(name_vec = c("sample", "library", "project", "medium", "percent.mito"))
  cl.categ  <- intersect(c(names(meta.tags), "Phase"), colnames(combined.obj@meta.data))
  (categ3Dplots <- unique(c(cl.categ, GetClusteringRuns(),  GetNamedClusteringRuns())))
  Plot3D.ListOfCategories(obj = combined.obj, cex = 1.75, ListOfCategories = categ3Dplots, annotate.by =  GetNamedClusteringRuns())
}


if (F) {
  TAP.genes.above.q1.25.for.3D.plots = c("ATF4", "BNIP3", "DDIT4", "EGR1", "ENO1", "ENO2", "FOS", "GAPDH", "HSPA8", "IER2", "JUN", "JUNB", "LDHA", "LDHB", "MARCKSL1", "PGAM1", "PGK1", "PKM", "SAT1", "TPI1", "PPP1R15A", "DDIT3", "SQSTM1", "P4HB")
  combined.obj <- gruffi::AddGOGeneList.manual(genes = TAP.genes.above.q1.25.for.3D.plots, GO = "TAP.genes.above.q1.25.for.3D.plots")
  Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = TAP.genes.above.q1.25.for.3D.plots, annotate.by =  GetNamedClusteringRuns()[1], cex = 2)
}


# End ------------------------
combined.obj <- RecallReduction(obj = combined.obj, reduction = "umap", dim = 2)
OutDir = getwd()


if (T) {
  Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = c('POLR2A', 'DCN', 'KAZN'), annotate.by =  GetNamedClusteringRuns()[1], cex = 2)
}


if (F) {
  combined.obj$cl.names.KnownMarkers.0.4
  all.genes$KAZN
  plot3D.umap.gene(obj = combined.obj, gene = 'N.RabV.N2c', AutoAnnotBy = GetNamedClusteringRuns()[1], dotsize = 2)
  plot3D.umap.gene(obj = combined.obj, gene = 'KAZN', AutoAnnotBy = GetNamedClusteringRuns()[1], dotsize = 2)
  plot3D.umap(obj = combined.obj, category = 'cl.names.KnownMarkers.0.4', AutoAnnotBy = 'cl.names.KnownMarkers.0.4', dotsize = 2)
  plot3D.umap(obj = combined.obj, category = 'cl.names.KnownMarkers.0.5', AutoAnnotBy = 'cl.names.KnownMarkers.0.5', dotsize = 2)

}

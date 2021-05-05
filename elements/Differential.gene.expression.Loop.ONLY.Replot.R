######################################################################
# Differential.gene.expression.Loop.ONLY.Replot.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.ONLY.Replot.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Differential.gene.expression.Loop.R")

# try(dev.off(), silent = T)

# Functions ------------------------
library(dplyr)
library(cowplot)

# Parameters ------------------------
# is.subclustering = p$"subclustering"
plotHeatmap <- TRUE

# Setup ------------------------
create_set_Original_OutDir()
create_set_SubDir("Differential.Gene.expression")

# if ( ww.variable.exists.and.true(p$"subclustering") ) {
#   create_set_SubDir("Differential.Gene.expression")
# } else {
#   create_set_OutDir(OutDirOrig,"Differential.Gene.expression")
# }
ParentDirDE <- OutDir


# Find ALL markers ------------------------------------------------------------------------
# p$'res.analyzed.DE' = c(.2,.5)

# plan("multiprocess", workers = 6)
# parallel.computing.by.future(workers_ = 6)

res.analyzed.DE = p$'res.analyzed.DE'
df.markers.all <- list.fromNames(res.analyzed.DE)
iprint("Resolutions analyzed: ", p$'res.analyzed.DE')
i = 2
for (i in 1:length(p$'res.analyzed.DE')) {

  res = p$'res.analyzed.DE'[i]
  create_set_OutDir(p0(ParentDirDE,ppp('res',res)))

  # Setup clustering identity for DE ------------------------------------
  p$'Ident.for.DEG' <-
    if (p$"cl.annotation"  == "ordered") {          GetOrderedClusteringRuns(res = res)
    # } else if (p$"cl.annotation"  == "named") {     GetNamedClusteringRuns(res = res)
    } else if (p$"cl.annotation"  == "simple") {    GetClusteringRuns(res = res)
    } else {print("not found")}

  stopifnot(p$'Ident.for.DEG' %in%  names(combined.obj@meta.data))
  Idents(combined.obj) <- p$'Ident.for.DEG'

  # Load DEG ------------------------------------

  df.markers <- combined.obj@misc$'df.markers'[[ppp('res',res)]]

  fname <- ppp('df.markers',res,'tsv')
  write.simple.tsv(df.markers, ManualName = fname)
  df.markers.all[[i]] <- df.markers

  clUMAP(ident = p$'Ident.for.DEG')
  PlotTopGenesPerCluster(obj = combined.obj, cl_res = res, nrGenes = p$'n.markers', order_by = p$"DEG.ranking"
                         , df_markers = combined.obj@misc$"df.markers"[[paste0("res.",res)]] )

  try(clUMAP(ident = ppp("cl.names.top.gene.res",res)))
  try(clUMAP(ident = ppp("cl.names.KnownMarkers",res)))


  plot.av.enrichment.hist = T
  if (plot.av.enrichment.hist) {
    df.markers.tbl <- as.tibble(df.markers)
    df.markers.tbl$cluster <- as.character(df.markers.tbl$cluster)
    p.deg.hist <- gghistogram(df.markers.tbl, x = "avg_log2FC",
                # add = "mean",
                rug = TRUE,
                color = "cluster", fill = "cluster",facet.by = 'cluster', xlim = c(0,3)
                # palette = c("#00AFBB", "#E7B800")
    ) + geom_vline(xintercept = 1)
    qqSave(p.deg.hist, w = 15, h = 15, title = ppp("Enrichment log2FC per cluster",res))

  }
} # for resolutions
write.simple.xlsx(named_list = df.markers.all )
p$"Cluster.Labels.Automatic" = F # so that it only runs 1x



# End ------------------------------------------------------------------------
create_set_Original_OutDir()


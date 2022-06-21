######################################################################
# Differential.gene.expression.Loop.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Differential.gene.expression.Loop.R")

# try(dev.off(), silent = T)

# Functions ------------------------
library(dplyr)
library(cowplot)
iprint('DefaultAssay',DefaultAssay(combined.obj))

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
i = 1

# p$'res.analyzed.DE' = c(0.3, 0.5)
for (i in 1:length(p$'res.analyzed.DE')) {
# p$"Cluster.Labels.Automatic" = T
# for (i in 1:1) {

  (res = p$'res.analyzed.DE'[i])
  create_set_OutDir(p0(ParentDirDE,ppp('res',res)))

  # Setup clustering identity for DE ------------------------------------
  p$'Ident.for.DEG' <-
    if (p$"cl.annotation"  == "ordered") {          GetOrderedClusteringRuns(res = res)
    # } else if (p$"cl.annotation"  == "named") {     GetNamedClusteringRuns(res = res)
    } else if (p$"cl.annotation"  == "simple") {    GetClusteringRuns(res = res)
    } else {print("not found")}

  stopifnot(p$'Ident.for.DEG' %in% names(combined.obj@meta.data))
  Idents(combined.obj) <- p$'Ident.for.DEG'

  # Find DEG ------------------------------------

  # Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE testing,
  # but could also miss features that are prefiltered

  tic(); df.markers <- FindAllMarkers(combined.obj, verbose = T
                                      , test.use = p$"test"
                                      , only.pos = p$"only.pos"
                                      , return.thresh	= p$"return.thresh"
                                      , min.pct = p$"min.pct"
                                      , min.diff.pct = p$"min.diff.pct"
                                      , min.cells.group = p$"min.cells.group"
                                      , min.cells.feature = p$"min.cells.feature"
                                      , logfc.threshold = p$"logfc.threshold"); toc();
  df.markers <- Add.DE.combined.score(df.markers) # , colLFC = 'avg_log2FC'

  combined.obj@misc$'df.markers'[[ppp('res',res)]] <- df.markers
  fname <- ppp('df.markers',res,'tsv')
  write.simple.tsv(df.markers, ManualName = fname)
  df.markers.all[[i]] <- df.markers

  clUMAP(ident = p$'Ident.for.DEG')
  PlotTopGenesPerCluster(obj = combined.obj, cl_res = res, nrGenes = p$'n.markers', order.by = p$"DEG.ranking"
                         , df_markers = combined.obj@misc$"df.markers"[[paste0("res.",res)]] )
  create_set_OutDir(ParentDirDE)

  if (p$"Cluster.Labels.Automatic") {

      # Auto Top genes ---------------------------------------------------------------
      combined.obj <- StoreAllMarkers(df_markers = df.markers, res = res)
      combined.obj <- AutoLabelTop.logFC(res = res);
      clUMAP(ident = ppp("cl.names.top.gene.res",res))

      # Hybrid ---------------------------------------------------------------
      combined.obj <- AddGOGeneList.manual(genes = p$"Cluster.Labels.Hybrid.Genes", GO = "BestMarkers")
      combined.obj <- AutoLabel.KnownMarkers(KnownMarkers = p$"Cluster.Labels.Hybrid.Genes", res = res)
      clUMAP(ident = ppp("cl.names.KnownMarkers",res))
  }

  plot.av.enrichment.hist = T
  if (plot.av.enrichment.hist) {
    df.markers.tbl <- as.tibble(df.markers)
    df.markers.tbl$cluster <- as.character(df.markers.tbl$cluster)
    p.deg.hist <- gghistogram(df.markers.tbl, x = "avg_log2FC",
                              title =  "Number of enriched genes per cluster",
                              subtitle =  "Binned by Log2(FC)",
                              caption =  paste(p$'Ident.for.DEG', "| v. line @ 2x FC"),
                              rug = TRUE,
                              color = "cluster", fill = "cluster",facet.by = 'cluster', xlim = c(0,3)
                              # palette = c("#00AFBB", "#E7B800")
    ) +
      geom_vline(xintercept = 1) +
      theme_linedraw()
    qqSave(p.deg.hist, w = 10, h = 6, title = ppp("Enrichment log2FC per cluster",res))

  }
} # for resolutions
write.simple.xlsx(named_list = df.markers.all )
p$"Cluster.Labels.Automatic" = F # so that it only runs 1x

iprint('DefaultAssay', DefaultAssay(combined.obj))

# End ------------------------------------------------------------------------
create_set_Original_OutDir()





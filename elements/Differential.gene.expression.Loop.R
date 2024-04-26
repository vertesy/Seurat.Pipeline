######################################################################
# Differential.gene.expression.Loop.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Differential.gene.expression.Loop.R")

# try(dev.off(), silent = T)
"https://github.com/vertesy/CON/issues/229"

"https://github.com/vertesy/CON/issues/229"

"https://github.com/vertesy/CON/issues/229"

"https://github.com/vertesy/CON/issues/229"


# Functions ------------------------------------------------
require(dplyr)
require(cowplot)
require(presto)
iprint('DefaultAssay',DefaultAssay(combined.obj))

# Parameters ------------------------------------------------
# is.subclustering = p$"subclustering"
plotHeatmap <- TRUE
rerun <- FALSE

# Setup ------------------------------------------------
create_set_Original_OutDir()
create_set_SubDir("Differential.Gene.expression.2")

ParentDirDE <- OutDir


# Checks ------------------------------------------------
stopifnot(exists('genes.ls'))


# Calculation: Find ALL markers ------------------------------------------------------------------------
res.analyzed.DE = p$'res.analyzed.DE'
df.markers.all <- list.fromNames(x = res.analyzed.DE)
iprint("Resolutions analyzed: ", p$'res.analyzed.DE')
i = 2

# p$'res.analyzed.DE' = c(0.3, 0.5)
for (i in 1:length(p$'res.analyzed.DE')) {
  (res = p$'res.analyzed.DE'[i])
  create_set_OutDir(p0(ParentDirDE,ppp('res', res)))

  # Setup clustering identity for DE ------------------------------------
  p$'Ident.for.DEG' <-
    if (p$"cl.annotation"  == "ordered") {          GetOrderedClusteringRuns(res = res)
    } else if (p$"cl.annotation"  == "simple") {    GetClusteringRuns(res = res)
    } else { print("not found") }

  stopifnot(p$'Ident.for.DEG' %in% names(combined.obj@meta.data))
  Idents(combined.obj) <- p$'Ident.for.DEG'

  # Find DEG ------------------------------------------------------------------------
  # Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE testing,
  # but could also miss features that are prefiltered

  if (rerun) {
    tic(); df.markers <- FindAllMarkers(combined.obj, verbose = T
                                        , test.use = p$"test"
                                        , only.pos = p$"only.pos"
                                        , return.thresh	= p$"return.thresh"
                                        , min.pct = p$"min.pct"
                                        , min.diff.pct = p$"min.diff.pct"
                                        , min.cells.group = p$"min.cells.group"
                                        # , min.cells.group = 50
                                        , logfc.threshold = p$"logfc.threshold"
                                        # , logfc.threshold = .5
                                        , max.cells.per.ident = p$"max.cells.per.ident"
                                        # , max.cells.per.ident = 100
    ); toc();

    df.markers <- Add.DE.combined.score(df.markers) # , colLFC = 'avg_log2FC'
    combined.obj@misc$'df.markers'[[ppp('res',res)]] <- df.markers
  } else {
    df.markers <- combined.obj@misc$'df.markers'[[ppp('res',res)]]
    if(is.null(df.markers)) { kollapse('df.markers is NULL. Entries in @misc$df.markers: ', names(combined.obj@misc$'df.markers')); stop()
    }
  }

  fname <- ppp('df.markers',res,'tsv')
  write.simple.tsv(df.markers, filename = fname)
  df.markers.all[[i]] <- df.markers
  xsave(df.markers, suffix = kpp("res", res))
} # for
create_set_OutDir(ParentDirDE)
xsave(combined.obj, suffix = kpp("w.DGEA", kpp('res', p$'res.analyzed.DE')))


# Stats and plots ------------------------------------------------------------------------
for (i in 1:length(p$'res.analyzed.DE')) {
  (res = p$'res.analyzed.DE'[i])
  create_set_OutDir(p0(ParentDirDE, ppp('res', res)))

  # Setup clustering identity for DE ------------------------------------
  p$'Ident.for.DEG' <-
    if (p$"cl.annotation"  == "ordered") {          GetOrderedClusteringRuns(res = res)
    } else if (p$"cl.annotation"  == "simple") {    GetClusteringRuns(res = res)
    } else { print("not found") }

  stopifnot(p$'Ident.for.DEG' %in% names(combined.obj@meta.data))
  Idents(combined.obj) <- p$'Ident.for.DEG'

  df.markers <- combined.obj@misc$'df.markers'[[ppp('res',res)]]
  if(is.null(df.markers)) { kollapse('df.markers is NULL. Entries in @misc$df.markers: ', names(combined.obj@misc$'df.markers')); stop() }

  clUMAP(ident = p$'Ident.for.DEG')
  # devtools::load_all(path = '~/GitHub/Packages/Seurat.utils');
  PlotTopGenesPerCluster(obj = combined.obj, cl_res = res, nrGenes = p$'n.markers', order.by = p$"DEG.ranking"
                         , df_markers = combined.obj@misc$"df.markers"[[paste0("res.",res)]])
  create_set_OutDir(ParentDirDE)

  if (p$"Cluster.Labels.Automatic") {

      # Auto Top genes ---------------------------------------------------------------
      combined.obj <- StoreAllMarkers(df_markers = df.markers, res = res)
      combined.obj <- AutoLabelTop.logFC(res = res);
      clUMAP(ident = ppp("cl.names.top.gene.res", res))

      # Hybrid ---------------------------------------------------------------
      combined.obj <- gruffi::AddGOGeneList.manual(genes = genes.ls$'Cluster.Labels.Hybrid.Genes.HiRes', GO = "BestMarkers")
      # combined.obj <- AutoLabel.KnownMarkers(KnownMarkers =  genes.ls$'Cluster.Labels.Hybrid.Genes.HiRes', res = res)
      # clUMAP(ident = ppp("cl.names.KnownMarkers",res))
  }

  plot.av.enrichment.hist = T
  if (plot.av.enrichment.hist) {
    df.markers.tbl <- as_tibble(df.markers)
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
    qqSave(ggobj = p.deg.hist, w = 10, h = 6, title = ppp("Enrichment log2FC per cluster",res))

  }
} # for resolutions
write.simple.xlsx(named_list = df.markers.all, rowname_column = 1 )
p$"Cluster.Labels.Automatic" = F # so that it only runs 1x

iprint('DefaultAssay', DefaultAssay(combined.obj))

# End ------------------------------------------------------------------------
create_set_Original_OutDir()





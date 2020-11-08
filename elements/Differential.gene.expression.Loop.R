######################################################################
# Differential.gene.expression.Loop.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.R')
# try(dev.off(), silent = T)

# Functions ------------------------
library(dplyr)

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
i=2
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

  # Find DEG ------------------------------------
  tic(); df.markers <- FindAllMarkers(combined.obj, verbose = T
                                      , test.use = p$"test"
                                      , only.pos = p$"only.pos"
                                      , return.thresh	= p$"return.thresh"
                                      , min.pct = p$"min.pct"
                                      , min.diff.pct = p$"min.diff.pct"
                                      , min.cells.group = p$"min.cells.group"
                                      , min.cells.feature = p$"min.cells.feature"
                                      , logfc.threshold = p$"logfc.threshold"); toc();
  df.markers <- Add.DE.combined.score(df.markers)

  combined.obj@misc$'df.markers'[[ppp('res',res)]] <- df.markers
  fname <- ppp('df.markers',res,'tsv')
  write.simple.tsv(df.markers, ManualName = fname)
  df.markers.all[[i]] <- df.markers

  clUMAP(ident = p$'Ident.for.DEG')
  PlotTopGenesPerCluster(obj = combined.obj, cl_res = res, nrGenes = p$'n.markers', order_by = p$"DEG.ranking"
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
} # for resolutions
write.simple.xlsx(named_list = df.markers.all )
p$"Cluster.Labels.Automatic" = F # so that it only runs 1x


# plotHeatmap ------------------------------------------------------------------------
if (plotHeatmap) {
  create_set_SubDir("heatmaps")
  p.hm <- list(
    assay = c('RNA', 'integrated')[2],
    slot = c('scale.data')
  )

  {
    NT <- c('GLRA2', 'GRM5', 'GRIN2B', 'GRIK2', 'GRID2', 'GRIA4', 'GRIA3', 'GRIA2', 'GRIA1', 'GABBR2', 'GABBR1', 'GABRG3', 'GABRG2', 'GABRG1', 'GABRB3', 'GABRB2', 'GABRA2', 'GABRA1', 'CNR1', 'HTR2C')
    p.hm <- list(assay = 'integrated', slot = 'scale.data', group.by = GetNamedClusteringRuns(res = p$'def_res'))
    hm.NT <- DoHeatmap(combined.obj, assay = p.hm$'assay', slot = p.hm$'slot', features = NT, group.by = p.hm$'group.by') # + NoLegend()
    ggsave2(plot = hm.NT, filename = ppp('hm.NT', flag.names_list.all.new(pl =  p.hm), 'pdf'))
  }

  p.hm$'rank'  <- c("combined.score", "avg_logFC", "p_val_adj", "p_val")[4]
  p.hm$'nrDEG' <- c(5, 10)[2]
  md.LogSettingsFromList(p.hm)
  {
    markers.use <- GetTopMarkers(dfDE = df.markers, n = p.hm$'nrDEG', order.by = p.hm$'rank')
    hm.NT <- DoHeatmap(combined.obj, assay = p.hm$'assay', slot = p.hm$'slot'
                       , features = as.character(markers.use), group.by = GetNamedClusteringRuns(res = p$'def_res')) # + NoLegend()
    ggsave2(plot = hm.NT, filename = ppp('hm.DEG', nrDEG,flag.names_list.all.new(pl= p.hm), 'pdf'), width = hA4, height = hA4)
  }
  create_set_Original_OutDir()
}


# End ------------------------------------------------------------------------
create_set_Original_OutDir()


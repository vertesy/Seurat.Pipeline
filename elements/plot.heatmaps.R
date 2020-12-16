######################################################################
# Plot.heatmaps.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.heatmaps.R')
# try(dev.off(), silent = T)

# Functions ------------------------
# library(dplyr)


# plotHeatmap ------------------------------------------------------------------------

create_set_SubDir("heatmaps")
p.hm <- list(
  assay = c('RNA', 'integrated')[1],
  group.by = GetNamedClusteringRuns(res = p$'def_res'),
  slot = c('scale.data')
)

# ------------------------------------------------------------------------
p.hm$'rank'  <- c("combined.score", "avg_logFC", "p_val_adj", "p_val")[1]
p.hm$'nrDEG' <- c(5, 10)[1]
md.LogSettingsFromList(p.hm)
{
  cl.markers <- unique(GetTopMarkers(dfDE = df.markers, n = p.hm$'nrDEG', order.by = p.hm$'rank'))
  hm.CM <- DoHeatmap(features = as.character(cl.markers), raster = FALSE
                     , object = combined.obj, assay = p.hm$'assay', slot = p.hm$'slot'
                     , group.by = GetNamedClusteringRuns(res = p$'def_res')) # + NoLegend()
  ggsave2(plot = hm.CM, filename = ppp('hm.DEG',  p.hm$'nrDEG',flag.names_list.all.new(pl = p.hm), 'pdf'), width = hA4, height = hA4)

  AvEx.CM <- AverageExpression(object = combined.obj, assays = 'RNA', slot = "data"
                               , features = cl.markers)
  pheatmap(mat = (AvEx.CM$RNA))
  wplot_save_this(plotname = "pheatmap.AvEx.CM", w = 10, h = 14)
  write.simple.tsv(iround(AvEx.CM$'RNA'))
}


# ------------------------------------------------------------------------
{
  ClassicMarkers
  "See for instance ('~/GitHub/Packages/Seurat.pipeline/Gene.Lists.Example.R')"

  hm.ClassicMarkers <- DoHeatmap(combined.obj, assay = p.hm$'assay', slot = p.hm$'slot', raster = FALSE
                                 , features = ClassicMarkers, group.by = p.hm$'group.by') # + NoLegend()
  ggsave2(plot = hm.ClassicMarkers, filename = ppp('hm.ClassicMarkers', flag.names_list.all.new(pl =  p.hm), 'pdf'))

  AvEx.ClassicMarkers <- AverageExpression(object = combined.obj, assays = 'RNA', slot = "data"
                                           , features = ClassicMarkers)
  pheatmap(mat = (AvEx.ClassicMarkers$RNA))
  wplot_save_this(plotname = "pheatmap.AvEx.ClassicMarkers", w = 10, h = 14)
  write.simple.tsv(iround(AvEx.ClassicMarkers$'RNA'))
}

create_set_Original_OutDir()


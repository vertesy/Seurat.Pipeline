######################################################################
# Differential.gene.expression.Loop.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/plot.heatmaps.R')
# try(dev.off(), silent = T)

# Functions ------------------------
library(dplyr)


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
  hm.CM <- DoHeatmap(features = as.character(cl.markers)
                     , object = combined.obj, assay = p.hm$'assay', slot = p.hm$'slot'
                     , group.by = GetNamedClusteringRuns(res = p$'def_res')) # + NoLegend()
  ggsave2(plot = hm.CM, filename = ppp('hm.DEG',  p.hm$'nrDEG',flag.names_list.all.new(pl = p.hm), 'pdf'), width = hA4, height = hA4)
  AvEx.CM <- AverageExpression(object = combined.obj, assays = 'RNA', slot = "data"
                               , features = cl.markers)
  pheatmap(mat = (AvEx.CM$RNA))
  wplot_save_this(plotname = "pheatmap.AvEx.CM", w = 10, h = 14)
  write.simple.tsv(iround(AvEx.CM$RNA))
}

# ------------------------------------------------------------------------
{
  genes.migration <- c("CXCR4", "EPHA3", "EPHA4", "EPHA5", "EPHB6", "KLF7", "PLXNA1",
                       "PLXNA4", "PLXNB1", "ROCK1", "SLIT1", "SLIT2", "PAK6", "SEMA3A",
                       "SEMA3F", "SEMA4F", "SEMA6C", "SEMA6D", "ROBO1", "ROBO2", "ROBO3",
                       "NRP2", "CNTN1", "DCC", "PPP3CA", "PPP3CB", "MAPK1", "MAPK3",
                       "DAB1", "ARX", "CNTNAP4", "LHX6", "SOX6", "ERBB4", "CXCL12",
                       "CXCL14", "NPY", "NXPH1", "NRG1", "NRG2", "NRG3", "NRG4")
  check.genes(list.of.genes = genes.migration, assay.slot = p.hm$'assay')

  hm.GM <- DoHeatmap(features = as.character(genes.migration)
                     , object = combined.obj, assay = p.hm$'assay', slot = p.hm$'slot'
                     , group.by = GetNamedClusteringRuns(res = p$'def_res')) # + NoLegend()
  ggsave2(plot = hm.GM, filename = ppp('hm.genes.migration', 'pdf'), width = hA4, height = hA4)
  AvEx.GM <- AverageExpression(object = combined.obj, assays = 'RNA', slot = "data"
                               , features = genes.migration)
  pheatmap(mat = (AvEx.GM$RNA))
  wplot_save_this(plotname = "pheatmap.AvEx.GM", w = 10, h = 14)
  write.simple.tsv(iround(AvEx.GM$RNA))
}

# ------------------------------------------------------------------------
{
  neuro.transmitters <- c('GLRA2', 'GRM5', 'GRIN2B', 'GRIK2', 'GRID2', 'GRIA4', 'GRIA3', 'GRIA2', 'GRIA1',
                          'GABBR2', 'GABBR1', 'GABRG3', 'GABRG2', 'GABRG1', 'GABRB3', 'GABRB2', 'GABRA2', 'GABRA1', 'CNR1', 'HTR2C')
  hm.NT <- DoHeatmap(combined.obj, assay = p.hm$'assay', slot = p.hm$'slot'
                     , features = neuro.transmitters, group.by = p.hm$'group.by') # + NoLegend()
  ggsave2(plot = hm.NT, filename = ppp('hm.NT', flag.names_list.all.new(pl =  p.hm), 'pdf'))
  AvEx.NT <- AverageExpression(object = combined.obj, assays = 'RNA', slot = "data"
                               , features = neuro.transmitters)
  pheatmap(mat = (AvEx.NT$RNA))
  wplot_save_this(plotname = "pheatmap.AvEx.NT", w = 10, h = 14)
  write.simple.tsv(iround(AvEx.NT$RNA))
}

create_set_Original_OutDir()


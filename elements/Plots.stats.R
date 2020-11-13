######################################################################
# Plot basic stats
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plots.stats.R')
try.dev.off()

# Setup ------------------------------------------------------------------------
create_set_Original_OutDir()
create_set_SubDir(ppp("01.Basic.Stats",p$"file.ext"))
stopifnot(exists('meta.tags'))

annot.clust <- GetClusteringRuns()
plnames.fixed.params = names(meta.tags) # , "RNA.model" , "organoid.fixed",  "organoid",

plots = c(annot.clust, plnames.fixed.params)
plot_list = list.fromNames(c(annot.clust, plnames.fixed.params))


# Make plots ------------------------------------------------------------------------
nrStat = length(annot.clust)
nrClass = length(plot_list)

for (pl in 1:length(plots)) {
  LB = (pl <=  nrStat || pl > nrClass)
  plot_list[[pl]] <- DimPlot(combined.obj, reduction = "umap", group.by = plots[pl], label = LB, repel = LB) +
    ggtitle(names(plot_list[pl]));
}
plot_list = plot_list[1:pl]

iprint('-------- ', length(plot_list), "plots are made")


# Draw and combine plots------------------------------------------------------------------------
lsplot.combinations <- iterBy.over(yourvec = plots, by = 4)
pl.names <-  c("clustering"
               , unlapply(lsplot.combinations, kpp)[-1])
names(lsplot.combinations) <- pl.names

for (i in 1:length(lsplot.combinations)) {
  plot.comination <- lsplot.combinations[[i]]
  fname <- ppp("UMAP", names(lsplot.combinations)[i], p$"file.ext")


  plotz = plot_grid(plotlist = plot_list[plot.comination], nrow = 2, ncol = 2, labels = LETTERS[1:l(plot.comination)]  )
  save_plot(filename = fname, plot = plotz, base_height = 12, ncol = 1, nrow = 1)
}
iprint('-------- ', length(plot_list), "plots are saved")

for (cr in GetClusteringRuns()) clUMAP(cr)

# Basic stats -----------------------------------
stats2plot <- intersect(p$"StatFeatures", colnames(combined.obj@meta.data))
stopifnot(l(stats2plot) > 0)

ggsave(FeaturePlot(combined.obj, min.cutoff = "q10", max.cutoff = "q90", reduction = 'umap', features = stats2plot)
       , filename = ppp("umap.StatMarkers", p$"file.ext"), width = hA4, height = wA4)

ggsave(FeaturePlot(combined.obj, min.cutoff = "q10", max.cutoff = "q90", reduction = 'tsne', features = stats2plot)
       , filename = ppp("tsne.StatMarkers", p$"file.ext"), width = hA4, height = wA4)


# Basic stat hexbinplots -----------------------------------
if (p$"plotHexBinStatPlots") {
  pl.hex <- foreach(i = 1:length(stats2plot)) %dopar% {
    plot_hexbin_meta(combined.obj, col = stats2plot[i], action = "median")
  }
  p.stat.hex = plot_grid(plotlist = pl.hex[1:4], nrow = 2, ncol = 2, labels = LETTERS[1:4]  )
  save_plot(filename = ppp("UMAPs.stat.hex", p$"file.ext"), plot = p.stat.hex, base_height = 12, ncol = 1, nrow = 1) #Figure 2

  pl.hex <- foreach(i = 1:length(plnames.fixed.params)) %dopar% {
    plot_hexbin_meta(combined.obj, col = plnames.fixed.params[i], action = "prop")
  }
  p.stat.hex2 = plot_grid(plotlist = pl.hex[1:4], nrow = 2, ncol = 2, labels = LETTERS[1:4]  )
  save_plot(filename = ppp("UMAPs.proportions.hex", p$"file.ext"), plot = p.stat.hex2, base_height = 12, ncol = 1, nrow = 1) #Figure 2
}



iprint('-------- StatMarker plots are saved')




#  -----------------------------------
#  -----------------------------------
if (p$"plotClusterPhylogeny") {
  # create_set_SubDir("ClusterPhylogeny")
  # Idents(combined.obj) <- GetNamedClusteringRuns(res = p$def_res);  # clUMAP()
  set.seed(p$seed)
  combined.obj <- BuildClusterTree(object = combined.obj, dims = 1:p$n.PC
                                   ,  assay = 'integrated', slot = 'scale.data',verbose = TRUE)
  PlotClusterTree(combined.obj)
  wplot_save_this("PlotClusterTree")

  if (FALSE) {
    # Plot on UMAP
    nr_cl <- max(as.numeric(Idents(combined.obj))) # Supercluster naming starts above
    x.ClusterTree <- Tool(object = combined.obj, slot = 'BuildClusterTree')
    nodeX <- nr_cl + 3
    plotsplit <- ColorDimSplit(combined.obj, node = nodeX) + ggtitle(label = ppp("Hierarchy node",nodeX))
  }
}
#  -----------------------------------

create_set_Original_OutDir()

######################################################################
# Principal Components PCA
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/PCA.R')
# try(dev.off(), silent = T)

# Functions ------------------------
# require(MarkdownReportsDev)

# Setup ------------------------
create_set_OutDir(OutDirOrig)
create_set_SubDir("PCA")


seu.plot.PC.var.explained(obj =  combined.obj)


# PCA.heatmap ------------------------
slotused <- if (p$'integrate.multiple') "integrated" else "RNA"
DefaultAssay(combined.obj) <- slotused

PCA.heatmap = T
if (PCA.heatmap) {
  AllPcs = list( c(1:12),c(13:24), c(25:30))
  for (pcDIMS in AllPcs) {
    iprint('Principal.Components', pcDIMS)
    PCs <- DimHeatmap(combined.obj, dims = pcDIMS, cells = 500,
                      raster = T, combine = F, fast = F)
    for (i in 1:length(PCs)) {PCs[[i]] <- PCs[[i]] + NoLegend()}

    PCs.c = plot_grid(plotlist = PCs, nrow = 4, ncol = 3,
                      labels = p0('PC ', pcDIMS), label_colour = "green", hjust = -2 )
    (fname = p0("Principal.Component.Loadings.",kppd(range(pcDIMS)),".pdf"))
    save_plot(plot = PCs.c, filename = fname, base_height = hA4, base_width = wA4)
    }
}

# PCA.plots ------------------------

PCA.plots = T
if (PCA.plots) {
  Idents(combined.obj) <- 'integrated_snn_res.0.3'
  plist.PCA = list.fromNames(2:13)
  for (PC in 2:13) {
    plist.PCA[[(PC - 1)]] <- (DimPlot(combined.obj, reduction = 'pca', label = T, dims = c(PC - 1, PC)) + NoLegend())
  }

  PCA.c = plot_grid(plotlist = plist.PCA, nrow = 4, ncol = 3,
                    labels = p0('PC ', 2:13) ) ;
  (fname = p0("Principal.Component.Dimensions.",kppd(range(2:13)),".pdf"))
  save_plot(plot = PCA.c, filename = fname, base_height = hA4, base_width = wA4)
}



# ------------------------
# ------------------------
# ------------------------
# ------------------------
# ------------------------


# End ------------------------------------------------------------------------
create_set_Original_OutDir()
DefaultAssay(combined.obj) <- "RNA"

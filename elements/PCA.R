######################################################################
# Principal Components PCA
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/PCA.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/PCA.R")
# try(dev.off(), silent = T)

# Functions ------------------------
# require(MarkdownReports)

# Setup ------------------------
create_set_OutDir(OutDirOrig)
create_set_SubDir("PCA")
axistextPCA = 6

seu.plot.PC.var.explained(obj =  combined.obj)


# PCA.heatmap ------------------------
slotused <- if (n.datasets > 1 & p$'integration' == 'CCA') "integrated" else "RNA"
DefaultAssay(combined.obj) <- slotused

PCA.heatmap = T
if (PCA.heatmap) {
  AllPcs = list( c(1:12),c(13:24), c(25:30))
  for (pcDIMS in AllPcs) {
    iprint('Principal.Components', pcDIMS)
    PCs <- DimHeatmap(combined.obj, dims = pcDIMS, cells = 500,
                      raster = T, combine = F, fast = F)
    for (i in 1:length(PCs)) {
      PCs[[i]] <- PCs[[i]] + NoLegend() +
      theme(axis.text.y = element_text(size = axistextPCA))
      }

    PCs.c = plot_grid(plotlist = PCs, nrow = 4, ncol = 3,
                      labels = p0('PC ', pcDIMS), label_colour = "green", hjust = -2 )
    (fname = p0("Principal.Component.Loadings.",kppd(range(pcDIMS)),".pdf"))
    save_plot(plot = PCs.c, filename = fname, base_height = hA4, base_width = wA4)
    }
}

# PCA.plots ------------------------

PCA.plots = T
if (PCA.plots) {
  Idents(combined.obj) <- identity.used[1]
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
if (T) {
  pcx <- combined.obj@reductions$pca@cell.embeddings
  idim(pcx)

  idx.num <-
    c( "nFeature_RNA", "percent.mito",
       "percent.ribo", "percent.AC.GenBank", "percent.AL.EMBL", "percent.LINC",
       "percent.MALAT1", "percent.antisense","log10.HGA_Markers",
       "S.Score",  "G2M.Score")
  idx.num <- intersect(idx.num, colnames(combined.obj@meta.data))
  cppx <- cor(cbind(pcx[,1:15], combined.obj@meta.data[, idx.num]))
  PCA.vs.Meta.cor.heatmap <- pheatmap::pheatmap(cppx, cutree_rows = 4, cutree_cols = 4)
  wplot_save_pheatmap(x = PCA.vs.Meta.cor.heatmap, width = 10 )

}

# ------------------------
# ------------------------
# ------------------------
# ------------------------


# End ------------------------------------------------------------------------
create_set_Original_OutDir()
DefaultAssay(combined.obj) <- "RNA"

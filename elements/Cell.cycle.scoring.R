######################################################################
# Cell Cycle Scoring
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Cell.cycle.scoring.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Cell.cycle.scoring.R")
# Based on https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
# try(dev.off(), silent = T)

# Parameters ------------------------
# plotCellFractionsBarplots = T
PlotRidgeplots = F
n.CC.genes = 16
plotCellFractionsBarplots = TRUE

stopifnot(all(names(meta.tags) %in% colnames(combined.obj@meta.data)))

# Setup ------------------------
create_set_OutDir(OutDirOrig, "Cell.cycle")
slotused <- if (n.datasets > 1 & p$'integration' == 'CCA') "integrated" else "RNA"

ccDir = "~/Dropbox (VBC)/Abel.IMBA/Metadata.D/gene.lists/cell_cycle_vignette_files/"
try(file.copy(from = ccDir, to = OutDir, recursive = T))
# regev_lab_cell_cycle_genes = read.simple.vec(ccDir,'regev_lab_cell_cycle_genes.txt')

# Metadata ------------------------
# regev_lab_cell_cycle_genes.extra = setdiff(regev_lab_cell_cycle_genes, unlist(cc.genes))
# multiFeaturePlot.A4(regev_lab_cell_cycle_genes.extra, obj = combined.obj)

# Parameters ------------------------
"Look at PCA for cell cycle effect in PCs"

combined.obj <- CellCycleScoring(combined.obj, s.features = cc.genes$'s.genes', g2m.features = cc.genes$'g2m.genes')
# view cell cycle scores and phase assignments

Idents(combined.obj) = "Phase"

# Read In ------------------------
(expr.q90.s.genes.found = names(sort(na.omit.strip(combined.obj@misc$'expr.q99'[ cc.genes$'s.genes']), decreasing = T) ) )
(expr.q90.g2m.genes.found = names(sort(na.omit.strip(combined.obj@misc$'expr.q99'[ cc.genes$'g2m.genes']), decreasing = T) ) )


expr.q90.s.genes.found = check.genes(expr.q90.s.genes.found, assay.slot = slotused)
expr.q90.g2m.genes.found = check.genes(expr.q90.g2m.genes.found, assay.slot = slotused)


Idents(combined.obj) <- "Phase"
Cell.cycle.umap.tsne.pca <- list(
  "umap" = DimPlot(combined.obj,  reduction = "umap"),
  "umap.r.0.3" = DimPlot(combined.obj,  reduction = "umap", group.by =  p$"res.MetaD.colname"),
  # "tsne" = DimPlot(combined.obj,  reduction = "tsne"),
  "pca.1.2" = DimPlot(combined.obj,  reduction = "pca", dims = c(1,2)),
  "pca.1.3" = DimPlot(combined.obj,  reduction = "pca", dims = c(1,3))
)
save4plotsA4(Cell.cycle.umap.tsne.pca)

try.dev.off()
CC.score = list(
  'G2M.Score' = qUMAP('G2M.Score'),
  'S.Score'   = qUMAP('S.Score')
)
save2plots.A4(CC.score)

clUMAP('Phase', label = F)

# CellFractionsBarplots ------------------------


if (plotCellFractionsBarplots) {
  pl.Fr = list(2)
  plotname = "Fraction.of.cell.cycle.stages.per.cluster"

  pl.Fr[[1]] <- scBarplot.CellFractions(obj = combined.obj, fill.by = "Phase", group.by = p$'res.MetaD.colname', downsample = F)
  pl.Fr[[2]] <- scBarplot.CellFractions(obj = combined.obj, group.by = "Phase", fill.by = "library", downsample = T)

  qqSaveGridA4(plotlist= pl.Fr, plots = 1:2, fname = ppp(plotname, "png"))
}


# CellFractionsBarplots ------------------------

if (plotCellFractionsBarplots) {
  pl.Fr = list.fromNames(names(meta.tags))

  # cols <- p$'res.MetaD.colname'
  cols <- GetNamedClusteringRuns(res = p$res.analyzed.DE[1])
  for (i in 1:l(meta.tags)) {
    pl.Fr[[i]] <- try(scBarplot.CellFractions(fill.by = names(meta.tags)[i], group.by = cols, downsample = T), silent = TRUE )
  }


  plotname = kpp("Fraction.of", names(meta.tags)[1], "and", names(meta.tags)[2], "per.cluster")
  qqSaveGridA4(plotlist= pl.Fr, plots = 1:2, fname = ppp(plotname, "pdf"))

  plotname = kpp("Fraction.of", names(meta.tags)[3], "and", names(meta.tags)[4], "per.cluster")
  qqSaveGridA4(plotlist= pl.Fr, plots = 3:4, fname = ppp(plotname, "pdf"))


}


# ------------------------
# ------------------------
# ------------------------

# PlotRidgeplots ------------------------
if (PlotRidgeplots) {
  iprint("ridgeplots ------------------------")
  rplot.s = RidgePlot(combined.obj, features = expr.q90.s.genes.found[1:6], ncol = 2)
  save_plot(plot = rplot.s, filename = "Cell.Cycle.S.1.pdf", base_height = hA4, base_width = wA4)

  rplot = RidgePlot(combined.obj, features = c("PCNA", "HELLS","BRIP1" ,"TOP2A","CENPF", "MKI67"), ncol = 2)
  save_plot(plot = rplot, filename = "Cell.Cycle.ridgeplot.pdf", base_height = hA4, base_width = wA4)

  # gm2 ------------------------
  rplot.g2m = RidgePlot(combined.obj, features =expr.q90.g2m.genes.found[1:6], ncol = 2)
  save_plot(plot = rplot.g2m, filename = "Cell.Cycle.g2m.1.pdf", base_height = hA4, base_width = wA4)

  # rplot.g2m = RidgePlot(combined.obj, features = expr.q90.g2m.genes.found[7:12], ncol = 2)
  # save_plot(plot = rplot.g2m, filename = "Cell.Cycle.g2m.2.pdf", base_height = hA4, base_width = wA4)
  #

}
# ------------------------







# Individual S and G2/M score plots ------------------------


if (FALSE) {

  # S ------------------------
  iprint("UMAP S ------------------------")
  multiFeaturePlot.A4(list.of.genes = expr.q90.s.genes.found[1:n.CC.genes], obj = combined.obj
                      , subdir = T, foldername = "S.genes")

  # g2m ------------------------
  iprint("UMAP g2m ------------------------")
  multiFeaturePlot.A4(list.of.genes = expr.q90.g2m.genes.found[1:n.CC.genes], obj = combined.obj
                      , subdir = T, foldername = "G2M.genes")

}



# End ------------------------------------------------------------------------
# Idents(combined.obj) = p$"res.MetaD.colname"
create_set_Original_OutDir()

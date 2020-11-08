######################################################################
# Cell Cycle Scoring
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Cell.cycle.scoring.R')
# Based on https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
# try(dev.off(), silent = T)

# Parameters ------------------------
PlotRidgeplots = F
n.CC.genes = 16

stopifnot(all(names(meta.tags) %in% colnames(combined.obj@meta.data)))

# Setup ------------------------
create_set_OutDir(OutDirOrig, "Cell.cycle")
slotused <- if (n.datasets > 1) "integrated" else "RNA"


ccDir = "~/Dropbox/Abel.IMBA/MetadataD/Gene.lists/cell_cycle_vignette_files/"
try(file.copy(from = ccDir, to = OutDir, recursive = T))
# regev_lab_cell_cycle_genes = read.simple.vec(ccDir,'regev_lab_cell_cycle_genes.txt')

# Metadata ------------------------
# regev_lab_cell_cycle_genes.extra = setdiff(regev_lab_cell_cycle_genes, unlist(cc.genes))
# multiFeaturePlot.A4(regev_lab_cell_cycle_genes.extra, obj = combined.obj)

# Parameters ------------------------
"Look at PCA for cell cycle effect in PCs"

combined.obj <- CellCycleScoring(combined.obj, s.features = cc.genes$'s.genes', g2m.features = cc.genes$'g2m.genes')
# view cell cycle scores and phase assignments
head(combined.obj)
Idents(combined.obj) = "Phase"

# Read In ------------------------
(expr.q90.s.genes.found = names(sort(na.omit.strip(combined.obj@misc$'expr.q90'[ cc.genes$'s.genes']), decreasing = T) ) )
(expr.q90.g2m.genes.found = names(sort(na.omit.strip(combined.obj@misc$'expr.q90'[ cc.genes$'g2m.genes']), decreasing = T) ) )


expr.q90.s.genes.found = check.genes(expr.q90.s.genes.found, assay.slot = slotused)
expr.q90.g2m.genes.found = check.genes(expr.q90.g2m.genes.found, assay.slot = slotused)


Idents(combined.obj) <- "Phase"
Cell.cycle.umap.tsne.pca <- list(
  "umap" = DimPlot(combined.obj,  reduction = "umap"),
  # "umap.r.0.3" = DimPlot(combined.obj,  reduction = "umap", group.by =  p$"res.MetaD.colname"),
  "tsne" = DimPlot(combined.obj,  reduction = "tsne"),
  "pca.1.2" = DimPlot(combined.obj,  reduction = "pca", dims = c(1,2)),
  "pca.1.3" = DimPlot(combined.obj,  reduction = "pca", dims = c(1,3))
)
save4umaps.A4(Cell.cycle.umap.tsne.pca)


CC.score = list(
  'G2M.Score' = qUMAP('G2M.Score'),
  'S.Score'   = qUMAP('S.Score')
)
save2umaps.A4(CC.score)

clUMAP('Phase')

# CellFractionsBarplots ------------------------
if (TRUE) {
  pl.Fr = list(2)
  plotname = "Fraction.of.cell.cycle.stages.per.cluster"

  pl.Fr[[1]] =CellFractionsBarplot2(obj = combined.obj, fill.by = "Phase", group.by = p$'res.MetaD.colname', downsample = F)
  pl.Fr[[2]] =CellFractionsBarplot2(obj = combined.obj, group.by = "Phase", fill.by = "sample", downsample = T)

  qqSaveGridA4(plotlist= pl.Fr, plots = 1:2, fname = ppp(plotname, "png"))
}



# CellFractionsBarplots ------------------------

if (plotCellFractionsBarplots) {
  pl.Fr = list(2)
  plotname = "Fraction.of.age.and.location.per.cluster"
  pl.Fr[[1]] = CellFractionsBarplot2(fill.by = "age", group.by = p$'res.MetaD.colname', downsample = T)
  pl.Fr[[2]] = CellFractionsBarplot2(fill.by = "location", group.by = p$'res.MetaD.colname', downsample = T)
  qqSaveGridA4(plotlist= pl.Fr, plots = 1:2, fname = ppp(plotname, "pdf"))

  plotname = "Fraction.of.samples.per.cluster"
  CellFractionsBarplot2(fill.by = "sample", group.by = p$'res.MetaD.colname')
  pl.Fr[[1]] = CellFractionsBarplot2(fill.by = "sample", group.by = p$'res.MetaD.colname', downsample = T)
  pl.Fr[[2]] = NULL
  qqSaveGridA4(plotlist= pl.Fr, plots = 1:2, fname = ppp(plotname, "pdf"))

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

expr.q90.s.genes.found
if (T) {

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
Idents(combined.obj) = p$"res.MetaD.colname"
create_set_Original_OutDir()

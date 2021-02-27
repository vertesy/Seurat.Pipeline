######################################################################
# Gene.expression.gene.lists.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R')


# create_set_SubDir("Gene.expression")
create_set_OutDir(OutDirOrig,"Gene.expression")
# OutDirOrig ="~/Dropbox/Abel.IMBA/AnalysisD/Bajaj/iN.migration..premRNA."


# Plot QC gene sets -----------------------------------

PlotTopGenes(obj = combined.obj)

# Highest.Expressed.Genes = names(head(sort(combined.obj@misc$expr.q90, decreasing = T), n = 16))
# multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, obj = combined.obj, intersectionAssay = "RNA", subdir =T)

if (l(p$Cluster.Labels.Hybrid.Genes)) {
  AutoNaming.Genes = gl$Cluster.Labels.Hybrid.Genes
  multiFeaturePlot.A4(list.of.genes = AutoNaming.Genes, obj = combined.obj, intersectionAssay = "RNA", subdir =T)
} else { iprint("p$Cluster.Labels.Hybrid.Genes IS NOT FOUND.")}


# ClassicMarkers.found <- check.genes(obj = combined.obj, list.of.genes = ClassicMarkers)
ClassicMarkers.plus = c(ClassicMarkers, LargeSubsetMarkers, QC.markers)
multiFeaturePlot.A4(list.of.genes = ClassicMarkers.plus, obj = combined.obj, subdir =T)
# multiFeaturePlot.A4(list.of.genes = ClassicMarkers, obj = combined.obj, plot.reduction = 'tsne', subdir =T)

AdditionalMarkers = c("BRN2", "TBR1", "TBR2", "VGLUT1", "VGLUT2", "NCAM", "L1CAM", "PSD95")
multiFeaturePlot.A4(list.of.genes = AdditionalMarkers, obj = combined.obj, subdir =T)


# End ------------------------------------------------------------------------
create_set_Original_OutDir()


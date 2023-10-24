######################################################################
# Gene.expression.gene.lists.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Gene.expression.gene.lists.R")

devtools::load_all(path = '~/GitHub/Packages/Seurat.utils');

# create_set_SubDir("Gene.expression")
create_set_OutDir(OutDirOrig,"Gene.expression")
# OutDirOrig ="~/Dropbox/Abel.IMBA/AnalysisD/Bajaj/iN.migration..premRNA."


# Plot QC gene sets -----------------------------------

PlotTopGenes(obj = combined.obj)

# Highest.Expressed.Genes = names(head(sort(combined.obj@misc$expr.q90, decreasing = T), n = 16))
# multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, obj = combined.obj, intersectionAssay = "RNA", subdir =T)

if (l(genes.ls$'Cluster.Labels.Hybrid.Genes')) {
  AutoNaming.Genes = genes.ls$'Cluster.Labels.Hybrid.Genes'
  multiFeaturePlot.A4(list.of.genes = AutoNaming.Genes, obj = combined.obj, intersectionAssay = "RNA", subdir =T)
} else { iprint("genes.ls$Cluster.Labels.Hybrid.Genes IS NOT FOUND.")}


OutDir
qUMAP(AutoNaming.Genes[1])
qUMAP(AutoNaming.Genes[2], )

# ClassicMarkers.found <- check.genes(obj = combined.obj, list.of.genes = ClassicMarkers)
ClassicMarkers.plus = c(genes.ls$'ClassicMarkers', genes.ls$'LargeSubsetMarkers', genes.ls$'QC.markers')
multiFeaturePlot.A4(list.of.genes = ClassicMarkers.plus, obj = combined.obj, subdir =T)
# multiFeaturePlot.A4(list.of.genes = ClassicMarkers, obj = combined.obj, plot.reduction = 'tsne', subdir =T)


AdditionalMarkers = c("BRN2", "TBR1", "TBR2", "VGLUT1", "VGLUT2", "NCAM", "L1CAM", "PSD95")
multiFeaturePlot.A4(list.of.genes = AdditionalMarkers, obj = combined.obj, subdir =T)



# Plot most variable genes from @var.features ------------------------------------------------------------------------
{

  Most.Variable.Genes <- head(combined.obj@assays$RNA@var.features, n=32)
  # toclip(combined.obj@assays$RNA@var.features)
  multiFeaturePlot.A4(list.of.genes = Most.Variable.Genes, obj = combined.obj, subdir =T)
}



if (p$'also.tSNE') multiFeaturePlot.A4(list.of.genes = ClassicMarkers.plus, obj = combined.obj, subdir =T, plot.reduction = 'tsne')


# End ------------------------------------------------------------------------
create_set_Original_OutDir()


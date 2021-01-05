######################################################################
# 00.Wrapper.for.sc.analysis.LOCAL.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/00.Wrapper.for.sc.analysis.LOCAL.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


# Functions ------------------------
source('~/GitHub/Packages/Seurat.pipeline/elements/Load.packages.local.R')
# sourceGitHub("Load.packages.CBE.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# Parameters ------------------------
source('~/GitHub/Packages/Seurat.pipeline/Parameters.example.R')
# sourceGitHub("Parameters.example.R" , repo = "Seurat.Pipeline", folder = ""); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/Gene.Lists.Example.R')
# sourceGitHub("Gene.Lists.Example.R" , repo = "Seurat.Pipeline", folder = ""); create_set_Original_OutDir()

# Setup ------------------------
OutDir = "~/sc.Analysis.Example/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "My.analysis.R")
OutDirOrig = OutDir

# Read In: Convert 10X output or load .RDS ------------------------
InputDir = "~/your/10X/folder"

MSeq.obj.fname = list.files(OutDir, pattern = "*.Rds")
MseqTabelExists = length(MSeq.obj.fname)

if (!MseqTabelExists) Convert10Xfolders(InputDir)
ls.Seurat <- LoadAllSeurats(InputDir, file.pattern =  "*.Rds")
(n.datasets = length(ls.Seurat))
(samples.short <- samples <- names(ls.Seurat))

ls.input <- ls.VarGenes <- list.fromNames(samples.short)

# Metadata ---------------------------
samples.short = "d80.AUT"
meta.tags <- list(
  'library' = samples.short,
  'project' = stringr::str_split_fixed(string = samples.short, pattern =  '\\.', n = 2)[,1],
  'sample' = stringr::str_split_fixed(string = samples.short, pattern =  '\\.', n = 3)[,2]
)


"The part below will not work out of the box."
for (i in 1:n.datasets ) {
  ls.Seurat[[i]] <- add.meta.fraction(col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-", obj = ls.Seurat[[i]])
  ls.Seurat[[i]] <- add.meta.fraction(col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS", obj = ls.Seurat[[i]])

  META = ls.Seurat[[i]]@meta.data
  ls.Seurat[[i]] <- AddMetaData(object = ls.Seurat[[i]], metadata = log10(META[, 'nCount_RNA']), col.name = 'log10.nCount_RNA')
  ls.Seurat[[i]] <- AddMetaData(object = ls.Seurat[[i]], metadata = log10(META[, 'nFeature_RNA']), col.name = 'log10.nFeature_RNA')
  ls.Seurat[[i]] <- AddMetaData(object = ls.Seurat[[i]], metadata = rep(samples.short[i], nCells), col.name = 'sample')
}

# Filtering & QC ------------------------
if (TRUE) PlotFilters(ls.obj = ls.Seurat); create_set_Original_OutDir()
source('~/GitHub/Packages/Seurat.pipeline/elements/Filtering.plots.3D.R'); create_set_Original_OutDir()
# sourceGitHub("Filtering.plots.3D.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

for (i in 1:n.datasets ) {
  iprint(names(ls.Seurat)[i], percentage_formatter(i/n.datasets, digitz = 2))
  sobj = ls.Seurat[[i]]
  sobj = subset(x = sobj, subset = `nFeature_RNA` > p$'thr.hp.nFeature_RNA' & `nFeature_RNA` < p$'thr.lp.nFeature_RNA')
  sobj = subset(x = sobj, subset = `percent.mito` > p$'thr.hp.mito' & `percent.mito` < p$'thr.lp.mito')
  sobj = subset(x = sobj, subset = `percent.ribo` > p$'thr.hp.ribo' & `percent.ribo` < p$'thr.lp.ribo')
  ls.Seurat[[i]] <- sobj
}; toc();

#
# How much of the gene sets overlap across objects? ------------------------
if (TRUE) source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.gene.set.overlaps.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.gene.set.overlaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R'); create_set_Original_OutDir()
# sourceGitHub("Gene.List.Overlap.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# Normalization ------------------------
for (i in 1:n.datasets ) {
  ls.Seurat[[i]] <- NormalizeData(object = ls.Seurat[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  ls.Seurat[[i]] <- FindVariableFeatures(object = ls.Seurat[[i]], mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR', nfeatures =10000)
}

source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.variable.Genes.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# If you want to downsample before integration ------------------------
if (FALSE) source('~/GitHub/Packages/Seurat.pipeline/elements/Dowsample.Seurat.Objects.R'); create_set_Original_OutDir()
# sourceGitHub("Dowsample.Seurat.Objects.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# Integrate ------------------------
useCCA = TRUE
if (useCCA) {
  tic(); anchors <- FindIntegrationAnchors(object.list = ls.Seurat, dims = 1:p$'n.CC', reduction = "rpca"); toc(); say();
  isave.RDS(anchors); isave.RDS(ls.Seurat, inOutDir = T)
  tic(); combined.obj <- IntegrateData(anchorset = anchors, dims = 1:p$'n.CC'); toc(); say()
} else {
 "use merge function"
}
DefaultAssay(combined.obj) <- "integrated"

## Perform integrated analysis --------------
combined.obj <- Calcq90Expression(obj = combined.obj)
# all.genes <- obj@misc$all.genes

# Scale / z-score norm ------------------------
for (assayX in c("RNA", "integrated")) {
  tic(); combined.obj <- ScaleData(combined.obj, assay = assayX, verbose = T, vars.to.regress = p$'variables.2.regress'); toc()
}


# t-SNE and Clustering ---------------------------------------------------------------------------
tic(); combined.obj <- RunPCA(combined.obj, npcs = p$'n.PC', verbose = T); toc()
tic(); combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
# tic(); combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions=3:2, reduction="umap"); toc()

tic(); combined.obj <- RunTSNE(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
tic(); combined.obj <- make_hexbin(combined.obj, nbins = 15, dimension_reduction = "UMAP"); toc()

tic(); combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
tic(); combined.obj <- FindClusters(combined.obj, resolution = p$'snn_res'); toc()
DefaultAssay(combined.obj) <- "RNA"
isave.RDS(combined.obj, inOutDir = T)


# Basic plots ------------------------
clUMAP()

qUMAP()

PlotTopGenes(obj = combined.obj)

source('~/GitHub/Packages/Seurat.pipeline/elements/Plots.stats.R'); create_set_Original_OutDir()
# sourceGitHub("Plots.stats.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R'); create_set_Original_OutDir()
# sourceGitHub("Gene.expression.gene.lists.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()



# Further calculations ------------------------

# Renumber clusters along a UMAP coordinate
source('~/GitHub/Packages/Seurat.pipeline/elements/Renumber.Clusters.R'); create_set_Original_OutDir()
# sourceGitHub("Renumber.Clusters.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# Calculate Differential gene expression
source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.R'); create_set_Original_OutDir()
# sourceGitHub("Differential.gene.expression.Loop.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()



# Further calculations and annotation ------------------------
source('~/GitHub/Packages/Seurat.pipeline/elements/PCA.R'); create_set_Original_OutDir()
# sourceGitHub("PCA.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/elements/Cell.cycle.scoring.R'); create_set_Original_OutDir()
# sourceGitHub("Cell.cycle.scoring.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()



# Plotting ------------------------
source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.3D.umaps.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.3D.umaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
combined.obj <- RecallReduction(obj = combined.obj, dim=2, reduction="umap")

source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.heatmaps.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.heatmaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# End ------------------------
if (TRUE) {  combined.obj@misc$p <- p; isave.RDS(object = combined.obj, inOutDir = T) } # Add a parameter list to a Seurat object's misc slot # seuSaveRds(object = combined.obj, use_Original_OutDir = T)
memory.biggest.objects()





# Not used, kept for reference ------------------------

if (FALSE) {
  source('~/GitHub/Packages/Seurat.pipeline/elements/JackStraw.plot.R'); create_set_Original_OutDir()
  # sourceGitHub("JackStraw.plot.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

  source('~/GitHub/Packages/Seurat.pipeline/elements/STRING.db.auto.cluster.annotation.R'); create_set_Original_OutDir()
  # sourceGitHub("STRING.db.auto.cluster.annotation.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
}
# ------------------------

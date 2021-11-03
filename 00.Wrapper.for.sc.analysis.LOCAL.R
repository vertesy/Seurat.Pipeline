######################################################################
# 00.Wrapper.for.sc.analysis.LOCAL.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/00.Wrapper.for.sc.analysis.LOCAL.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


"This is not a compete dummy pipeline, just steps of Seurat.Pipeline in order, many intermediate steps are missing"
"This is not a compete dummy pipeline, just steps of Seurat.Pipeline in order, many intermediate steps are missing"
"This is not a compete dummy pipeline, just steps of Seurat.Pipeline in order, many intermediate steps are missing"

# Functions ------------------------
# sourceGitHub is a function in CodeAndRoll.R
source('~/GitHub/Packages/Seurat.pipeline/elements/Load.packages.local.R')
# sourceGitHub("Load.packages.CBE.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "My.analysis.R")
OutDirOrig = OutDir


# Parameters ------------------------
source('~/GitHub/Packages/Seurat.pipeline/Parameters.example.R')
# sourceGitHub("Parameters.example.R" , repo = "Seurat.Pipeline", folder = ""); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/Gene.Lists.Example.R')
# sourceGitHub("Gene.Lists.Example.R" , repo = "Seurat.Pipeline", folder = ""); create_set_Original_OutDir()

# Read In ------------------------
if (!MseqTabelExists) Convert10Xfolders(InputDir)
ls.Seurat <- LoadAllSeurats(InputDir, file.pattern =  "*.Rds")
(samples <- names(ls.Seurat))
ls.input <- ls.VarGenes <- list.fromNames(samples)
(n.datasets = length(ls.Seurat))


# meta.tags ------------------------
# Adjust this to each dataset
meta.tags <- list(
  'library' = samples,
  'project' = stringr::str_split_fixed(string = samples, pattern =  '\\.', n = 2)[,1],
  'sample' = stringr::str_split_fixed(string = samples, pattern =  '\\.', n = 3)[,2]
)
lapply(meta.tags, unique)

# QC ------------------------

# Filtering ------------------------
PlotFilters(ls.obj = ls.Seurat)
source('~/GitHub/Packages/Seurat.pipeline/elements/Filtering.plots.3D.R'); create_set_Original_OutDir()
# sourceGitHub("Filtering.plots.3D.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# Filter ------------------------
for (i in 1:n.datasets ) {
  iprint(names(ls.Seurat)[i], percentage_formatter(i/n.datasets, digitz = 2))
  sobj = ls.Seurat[[i]]
  sobj = subset(x = sobj, subset = `nFeature_RNA` > p$'thr.hp.nFeature_RNA' & `nFeature_RNA` < p$'thr.lp.nFeature_RNA')
  sobj = subset(x = sobj, subset = `percent.mito` > p$'thr.hp.mito' & `percent.mito` < p$'thr.lp.mito')
  sobj = subset(x = sobj, subset = `percent.ribo` > p$'thr.hp.ribo' & `percent.ribo` < p$'thr.lp.ribo')
  ls.Seurat[[i]] <- sobj
}; toc();


# Gene set comparisons ------------------------
source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.gene.set.overlaps.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.gene.set.overlaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R'); create_set_Original_OutDir()
# sourceGitHub("Gene.List.Overlap.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# Downsample individual datasets if needed ------------------------
if (FALSE) {
  source('~/GitHub/Packages/Seurat.pipeline/elements/Dowsample.Seurat.Objects.R'); create_set_Original_OutDir()
  # sourceGitHub("Dowsample.Seurat.Objects.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
}

# Integrate ------------------------
tic(); anchors <- FindIntegrationAnchors(object.list = ls.Seurat, dims = 1:p$'n.CC'); toc(); say();
isave.RDS(anchors); isave.RDS(ls.Seurat)
tic(); combined.obj <- IntegrateData(anchorset = anchors, dims = 1:p$'n.CC', features.to.integrate = features.2.integrate); toc(); say()

# Save meta.data ------------------------
combined.obj <- calc.q90.Expression.and.set.all.genes(obj = combined.obj)  # Sets combined.obj@misc$'all.genes' and 'all.genes' variables.
combined.obj@misc$'n.datasets'  <- n.datasets
combined.obj@misc$'meta.tags'   <- meta.tags
combined.obj@misc$'p'           <- p

isave.RDS(combined.obj, inOutDir = T)

# Clustering & co ------------------------
DefaultAssay(combined.obj) <- "integrated"
tic(); combined.obj <- ScaleData(combined.obj, verbose = T, vars.to.regress = p$'variables.2.regress'); toc()
tic(); combined.obj <- RunPCA(combined.obj, npcs = p$'n.PC', verbose = T); toc()
tic(); combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions = 3:2, reduction = "umap"); toc()
tic(); combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
tic(); combined.obj <- FindClusters(combined.obj, resolution = p$'snn_res'); toc()
isave.RDS(combined.obj, inOutDir = T)
tic(); combined.obj <- RunTSNE(combined.obj, reduction = "pca", dims = 1:p$'n.PC', method = "FIt-SNE"); toc() # https://www.rdocumentation.org/packages/Seurat/versions/4.0.1/topics/RunTSNE



# Basic plots ------------------------
source('~/GitHub/Packages/Seurat.pipeline/elements/Plots.stats.R'); create_set_Original_OutDir()
# sourceGitHub("Plots.stats.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R'); create_set_Original_OutDir()
# sourceGitHub("Gene.expression.gene.lists.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# Further calculations ------------------------

# Renumber clusters along a UMAP coordinate
# 'You need to manually redefine p$"Reorder.Dim".  -1 means right to left aling umap dimension 1.'
p$"Reorder.Dim" <- 1
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

source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.heatmaps.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.heatmaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R'); create_set_Original_OutDir()
# sourceGitHub("Plot.variable.Genes.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# End ------------------------
if (TRUE) { isave.RDS(object = combined.obj, inOutDir = T) } # Add a parameter list to a Seurat object's misc slot # seuSaveRds(object = combined.obj, use_Original_OutDir = T)
memory.biggest.objects()


# ContinuePreviousAnalysis ------------------------
if (FALSE) {
  rdsDir = "~/path/to/_RDS/"
  combined.obj <- read_rds(kpps(rdsDir, 'combined.obj__2021.04.22_14.10.Rds.gz'))

  recall.parameters()
  recall.all.genes()
  recall.meta.tags.n.datasets()
}

# Not used, kept for reference ------------------------

if (FALSE) {
  source('~/GitHub/Packages/Seurat.pipeline/elements/JackStraw.plot.R'); create_set_Original_OutDir()
  # sourceGitHub("JackStraw.plot.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

  source('~/GitHub/Packages/Seurat.pipeline/elements/STRING.db.auto.cluster.annotation.R'); create_set_Original_OutDir()
  # sourceGitHub("STRING.db.auto.cluster.annotation.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
}
# ------------------------

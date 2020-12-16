######################################################################
# 00.Wrapper.for.sc.analysis.WEB.R
######################################################################
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/00.Wrapper.for.sc.analysis.LOCAL.R')
# source('~/GitHub/Packages/Seurat.pipeline/00.Wrapper.for.sc.analysis.WEB.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


"This is not a compete dummy pipeline, just steps of Seurat.Pipeline in order, many intermediate steps are missing"
"This is not a compete dummy pipeline, just steps of Seurat.Pipeline in order, many intermediate steps are missing"
"This is not a compete dummy pipeline, just steps of Seurat.Pipeline in order, many intermediate steps are missing"

# Functions ------------------------
# sourceGitHub is a function in CodeAndRoll.R
sourceGitHub("Load.packages.CBE.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Load.packages.local.R')


# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "My.analysis.R")
OutDirOrig = OutDir

# Metadata ------------------------

# Parameters ------------------------
sourceGitHub("Parameters.example.R" , repo = "Seurat.Pipeline", folder = ""); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/Parameters.example.R')

sourceGitHub("Gene.Lists.Example.R" , repo = "Seurat.Pipeline", folder = ""); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/Gene.Lists.Example.R')

# Read In ------------------------

# QC ------------------------

# Filtering ------------------------

sourceGitHub("Filtering.plots.3D.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Filtering.plots.3D.R'); create_set_Original_OutDir()

# ------------------------
sourceGitHub("Plot.gene.set.overlaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.gene.set.overlaps.R'); create_set_Original_OutDir()

sourceGitHub("Gene.List.Overlap.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R'); create_set_Original_OutDir()

sourceGitHub("Dowsample.Seurat.Objects.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Dowsample.Seurat.Objects.R'); create_set_Original_OutDir()


# Basic plots ------------------------
sourceGitHub("Plots.stats.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plots.stats.R'); create_set_Original_OutDir()

sourceGitHub("Gene.expression.gene.lists.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R'); create_set_Original_OutDir()


# Further calculations ------------------------

# Renumber clusters along a UMAP coordinate
sourceGitHub("Renumber.Clusters.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Renumber.Clusters.R'); create_set_Original_OutDir()

# Calculate Differential gene expression
sourceGitHub("Differential.gene.expression.Loop.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Differential.gene.expression.Loop.R'); create_set_Original_OutDir()



# Further calculations and annotation ------------------------
sourceGitHub("PCA.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/PCA.R'); create_set_Original_OutDir()

sourceGitHub("Cell.cycle.scoring.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Cell.cycle.scoring.R'); create_set_Original_OutDir()



# Plotting ------------------------
sourceGitHub("Plot.3D.umaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.3D.umaps.R'); create_set_Original_OutDir()

sourceGitHub("Plot.heatmaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.heatmaps.R'); create_set_Original_OutDir()

sourceGitHub("Plot.variable.Genes.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R'); create_set_Original_OutDir()


# End ------------------------





# Not used, kept for reference ------------------------

if (FALSE) {
  sourceGitHub("JackStraw.plot.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
  # source('~/GitHub/Packages/Seurat.pipeline/elements/JackStraw.plot.R'); create_set_Original_OutDir()

  sourceGitHub("STRING.db.auto.cluster.annotation.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
  # source('~/GitHub/Packages/Seurat.pipeline/elements/STRING.db.auto.cluster.annotation.R'); create_set_Original_OutDir()
}
# ------------------------

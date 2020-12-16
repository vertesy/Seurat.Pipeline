######################################################################
# 00.Wrapper.for.sc.analysis.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/00.Wrapper.for.sc.analysis.R)
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/00.Wrapper.for.sc.analysis.R')

# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
sourceGitHub("Load.packages.CBE.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# sourceGitHub("Load.packages.local.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# source('~/GitHub/Packages/Seurat.pipeline/elements/Load.packages.local.R')


# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "")
OutDirOrig = OutDir

# Metadata ------------------------

# Parameters ------------------------


# Read In ------------------------

# QC ------------------------

# ------------------------
sourceGitHub("Filtering.plots.3D.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# sourceGitHub("Filtering.plots.3D.multiplex.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
# sourceGitHub("Filtering.plots.3D.solo.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# ------------------------
sourceGitHub("plot.gene.set.overlaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("Gene.List.Overlap.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("Dowsample.Seurat.Objects.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# ------------------------
sourceGitHub("Plots.stats.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("Gene.expression.gene.lists.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()


# ------------------------
sourceGitHub("Renumber.Clusters.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("Differential.gene.expression.Loop.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# ------------------------

sourceGitHub("PCA.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("Cell.cycle.scoring.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

sourceGitHub("plot.3D.umaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("plot.heatmaps.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
sourceGitHub("Plot.variable.Genes.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

# ------------------------

# Not used ------------------------

if (FALSE) {
  sourceGitHub("JackStraw.plot.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()
  sourceGitHub("STRING.db.auto.cluster.annotation.R" , repo = "Seurat.Pipeline"); create_set_Original_OutDir()

}


# ------------------------


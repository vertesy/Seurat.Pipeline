######################################################################
# Pairwise.Std.var.correlation.R
######################################################################
# source ('~/GitHub/Packages/Seurat.pipeline/elements/Pairwise.Std.var.correlation.R')

# Setup ------------------------------------------------------------------------
ls.VarGenes.top20 = list.fromNames(names(ls.Seurat))
create_set_OutDir(OutDirOrig, "variable.genes")

# Plot ------------------------------------------------------------------------



tic(); for (i in 1:n.datasets ) { print(i)
  ls.VarGenes.top20[[i]] <- head(VariableFeatures(ls.Seurat[[i]]), 20) # Identify the 10 most highly variable genes
  plot1 <- VariableFeaturePlot(ls.Seurat[[i]])
  ggsave(LabelPoints(
    plot = plot1, points = ls.VarGenes.top20[[i]], repel = TRUE)
    , filename = ppp("Var.genes",samples[i],'pdf') )
}; toc()
# wvenn(ls.VarGenes.top20)

pairwise.scatters=T
if (pairwise.scatters && !p$'use.SCTransform') {
  topN =100
  ls.variance.standardized =
    lapply(
      lapply(
        lapply(
          lapply(
            lapply(ls.Seurat,
                   FUN = HVFInfo),
            FUN = `[`, i='variance.standardized'),
          FUN = as.named.vector, WhichDimNames = 1),
        FUN = sort, decreasing = TRUE),
      FUN = head, n = topN)
  names(ls.variance.standardized) = samples

  plotname=ppp("Pairwise Correlation of Standardized Variance\ntopN", topN)
  pairs(ls.variance.standardized, main=plotname, upper.panel = panel.cor.pearson)
  # pairs(ls.variance.standardized, main=plotname ,
  #       col=rgb(0,0,0,0.25), pch=18) # , upper.panel = panel.cor.pearson
  #
  wplot_save_this(plotname, w = 10)
}

# End ------------------------------------------------------------------------

create_set_Original_OutDir()

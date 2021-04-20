######################################################################
# Plot.variable.Genes.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Plot.variable.Genes.R")

# Setup ------------------------------------------------------------------------
ls.VarGenes.top20 = list.fromNames(names(ls.Seurat))
create_set_OutDir(OutDirOrig, "variable.genes")

pairwise.scatters=T

# Plot ------------------------------------------------------------------------
tic(); for (i in 1:n.datasets ) { print(i)
  ls.VarGenes.top20[[i]] <- head(VariableFeatures(ls.Seurat[[i]]), 20) # Identify the 10 most highly variable genes
  plot1 <- VariableFeaturePlot(ls.Seurat[[i]])
  ggsave(LabelPoints(
    plot = plot1, points = ls.VarGenes.top20[[i]], repel = TRUE)
    , filename = ppp("Var.genes",samples[i],'pdf') )
}; toc()
# wvenn(ls.VarGenes.top20)


# Plot pairwise.scatters------------------------------------------------------------------------
if (pairwise.scatters ) {
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

  plotname=ppp("Pairwise Correlation of Standardized Variance of top", topN,"genes.")
  size <- round(length(ls.Seurat)/2)+5
  try.dev.off()
  pdf(file = "Genes standardized variance correlation across datasets.pdf",width = size, height = size)
  x <- pairs(ls.variance.standardized, main=plotname, upper.panel = panel.cor.pearson)
  try.dev.off()

}

# Heatmap ------------------------------------------------------------------------
# if (FALSE) {
#   ls.variance.standardized.1000 =
#     lapply(
#       lapply(
#         lapply(
#           lapply(
#             lapply(ls.Seurat,
#                    FUN = HVFInfo),
#             FUN = `[`, i='variance.standardized'),
#           FUN = as.named.vector, WhichDimNames = 1),
#         FUN = sort, decreasing = TRUE),
#       FUN = head, n = topN)
#   names(ls.variance.standardized.1000) = samples
#   df.x <- list2fullDF.byNames(ls.variance.standardized.1000)
#
#
#   pho <- pheatmap::pheatmap(df.x, show_rownames = F)
#
# }

# End ------------------------------------------------------------------------

create_set_Original_OutDir()

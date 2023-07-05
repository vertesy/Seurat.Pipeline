######################################################################
# Plot.variable.Genes.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Plot.variable.Genes.R")

# Setup ------------------------------------------------------------------------
ls.VarGenes.top20 = list.fromNames(names(ls.Seurat))
create_set_OutDir(OutDirOrig, "variable.genes")

pairwise.scatters=T

# n.datasets = 1
# samples = "sc4"


if(!exists('n.datasets')) n.datasets <- length(ls.Seurat)
if(!exists('samples')) samples <- paste0("Sample.", 1:length(ls.Seurat))


# Plot ------------------------------------------------------------------------
tic();
# Identify the 10 most highly variable genes in each dataset
ls.VarGenes.top20 <- lapply(ls.Seurat, function(x) head(VariableFeatures(x), 20))

# Create a variable feature plot for each dataset, labeling the top 20 most variable genes
for (i in 1:n.datasets) {
  plot1 <- VariableFeaturePlot(ls.Seurat[[i]])
  ggsave(LabelPoints(
    plot = plot1, points = ls.VarGenes.top20[[i]], repel = TRUE)
    , width = 7, h = 5
    , filename = ppp("Var.genes",samples[i],'pdf') )
}; toc()




# Plot pairwise.scatters------------------------------------------------------------------------
# if (pairwise.scatters && n.datasets >1) {
#   topN =100
#   ls.variance.standardized =
#     lapply(
#       lapply(
#         lapply(
#           lapply(
#             lapply(ls.Seurat,
#                    FUN = HVFInfo),
#             FUN = `[`, i='variance.standardized'),
#           FUN = as.named.vector.df, WhichDimNames = 1),
#         FUN = sort, decreasing = TRUE),
#       FUN = head, n = topN)
#   names(ls.variance.standardized) = samples
#
#   plotname=ppp("Pairwise Correlation of Standardized Variance of top", topN,"genes.")
#   size <- round(length(ls.Seurat)/2)+5
#   try.dev.off()
#   pdf(file = "Genes standardized variance correlation across datasets.pdf",width = size, height = size)
#   x <- pairs(ls.variance.standardized, main=plotname, upper.panel = panel.cor.pearson)
#   try.dev.off()
# }

if (pairwise.scatters && n.datasets > 1) {
  topN = 100
  ls.variance.standardized <- purrr::map(ls.Seurat, ~ {
    genes <- rownames(HVFInfo(.)$variance.standardized)
    data <- cbind(genes, unlist(HVFInfo(.)$variance.standardized))
    as.named.vector.df(data, WhichDimNames = 1)
  }) %>%
    purrr::map(~ sort(., decreasing = TRUE)) %>%
    purrr::map(~ head(., topN))
  names(ls.variance.standardized) <- samples

  plotname = ppp("Pairwise Correlation of Standardized Variance of top", topN, "genes.")
  size <- round(length(ls.Seurat) / 2) + 5
  try.dev.off()
  pdf(file = "Genes standardized variance correlation across datasets.pdf",width = size, height = size)
  x <- pairs(ls.variance.standardized, main=plotname, upper.panel = panel.cor.pearson)
  try.dev.off()
}

# Heatmap ------------------------------------------------------------------------

# End ------------------------------------------------------------------------

create_set_Original_OutDir()

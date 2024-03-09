######################################################################
# Plot.variable.Genes.Combined.R
######################################################################
file.remove('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.Combined.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Plot.variable.Genes.Combined.R")


VariableFeaturePlot
# Functions ------------------------------------------------------------------------


# Setup ------------------------------------------------------------------------
create_set_OutDir(OutDirOrig, "variable.genes")

pairwise.scatters = T




# Plot ------------------------------------------------------------------------
# Identify the 10 most highly variable genes in each dataset
ls.VarGenes.top20 <- lapply(ls.Seurat, function(x) head(VariableFeatures(x), 20))

# Create a variable feature plot for each dataset, labeling the top 20 most variable genes
for (i in 1:length(ls.Seurat)) {
  plot1 <- Seurat::VariableFeaturePlot(ls.Seurat[[i]])
  ggplot2::ggsave(LabelPoints(
    plot = plot1, points = ls.VarGenes.top20[[i]], repel = TRUE)
    , width = 7, h = 5
    , filename = ppp("Var.genes",samples[i],'pdf') )
};




# Plot pairwise.scatters------------------------------------------------------------------------
if (pairwise.scatters && length(ls.Seurat) > 1) {
  topN = 100
  ls.variance.standardized <- purrr::map(ls.Seurat, ~ {
    genes <- rownames(HVFInfo(.)$variance.standardized)
    data <- cbind(genes, unlist(HVFInfo(.)$variance.standardized))
    as.named.vector.df(data, WhichDimNames = 1)
  }) |>
    purrr::map(~ sort(., decreasing = TRUE)) |>
    purrr::map(~ head(., topN))
  names(ls.variance.standardized) <- samples

  plotname = ppp("Pairwise Correlation of Standardized Variance of top", topN, "genes.")
  size <- round(length(ls.Seurat) / 2) + 5
  try.dev.off()
  grDevices::pdf(file = "Genes standardized variance correlation across datasets.pdf",width = size, height = size)
  x <- graphics::pairs(ls.variance.standardized, main=plotname, upper.panel = panelCorPearson)
  try.dev.off()
}

# Heatmap ------------------------------------------------------------------------

# End ------------------------------------------------------------------------

create_set_Original_OutDir()


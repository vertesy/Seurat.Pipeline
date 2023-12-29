######################################################################
# Plot.variable.Genes.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Plot.variable.Genes.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Plot.variable.Genes.R")



# Functions ------------------------------------------------------------------------
# from rocinante
panel.cor.pearson <- function(x, y, digits = 2, prefix = "", cex.cor = 2, method = "pearson") { # A function to display correlation values for pairs() function. Default is pearson correlation, that can be set to  "kendall" or "spearman".
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = method, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  test <- cor.test(x, y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  
  text(0.5, 0.5, txt, cex = cex * r)
  text(.8, .8, Signif, cex = cex,  col = 2)
}


# Setup ------------------------------------------------------------------------
ls.VarGenes.top20 = list.fromNames(names(ls.Seurat))
create_set_OutDir(OutDirOrig, "variable.genes")

pairwise.scatters = T


if(!exists('samples')) samples <- paste0("Sample.", 1:length(ls.Seurat))


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
  x <- graphics::pairs(ls.variance.standardized, main=plotname, upper.panel = panel.cor.pearson)
  try.dev.off()
}

# Heatmap ------------------------------------------------------------------------

# End ------------------------------------------------------------------------

create_set_Original_OutDir()


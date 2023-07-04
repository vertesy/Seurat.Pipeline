######################################################################
# Gene.List.Overlap.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Gene.List.Overlap.R")

# Functions ------------------------
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))


# Read In ------------------------
ls_genes  <- lapply(ls.Seurat, rownames);

# Calculate Jaccard ------------------------
PairwiseJaccardIndices <- jPairwiseJaccardIndexList(lsG = ls_genes)

# ------------------------
colnames(PairwiseJaccardIndices)

ph.Jac.wo <- pheatmap::pheatmap(PairwiseJaccardIndices, display_numbers = T, cutree_rows = 3, cutree_cols = 3)
wplot_save_pheatmap(ph.Jac.wo, width = max(round(n.datasets*0.3), 8))

nrGenesPerDataset <- unlapply(ls_genes, length)
qbarplot(nrGenesPerDataset
         , subtitle = "Not actually filtered on this thr.", hline = 20000
         , xlab.angle = 45)

# End ------------------------
create_set_Original_OutDir()



######################################################################
# Gene.List.Overlap.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Gene.List.Overlap.R")

# Functions ------------------------
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))

n_cuts <- min(4, length(ls.Seurat))-1

# Read In ------------------------
ls_genes  <- lapply(ls.Seurat, rownames);

# Calculate Jaccard ------------------------
PairwiseJaccardIndices <- jPairwiseJaccardIndexList(lsG = ls_genes)

# ------------------------
colnames(PairwiseJaccardIndices)


ph.Jac.wo <- pheatmap::pheatmap(PairwiseJaccardIndices, display_numbers = T, cutree_rows = , cutree_cols = n_cuts)
wplot_save_pheatmap(ph.Jac.wo, width = max(round(length(ls.Seurat)*0.3), 8))

nrGenesPerDataset <- unlapply(ls_genes, length)
qbarplot(nrGenesPerDataset, label = nrGenesPerDataset
         , subtitle = "Not actually filtered on this thr.", hline = 20000
         , xlab.angle = 45)

# End ------------------------
create_set_Original_OutDir()



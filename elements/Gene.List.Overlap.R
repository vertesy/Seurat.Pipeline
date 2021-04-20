######################################################################
# Gene.List.Overlap.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R')
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Gene.List.Overlap.R")
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))

# source('~/Github/TheCorvinas/R/DatabaseLinke.r')
require(MarkdownReportsDev)

# Setup ------------------------
# OutDir = "/Volumes/abel/cbehome/Dropbox/Abel.IMBA/AnalysisD/Abel/SEO2/INTEGR~1/Gene.List.Overlap"
# setup_MarkdownReports(OutDir = OutDir, scriptname = "Gene.List.Overlap.R")
# OutDirOrig = OutDir

# Metadata ------------------------
# Parameters ------------------------
# Read In ------------------------
# ls_genes <- read_rds(file = "/Volumes/abel/cbehome/Dropbox/Abel.IMBA/AnalysisD/Abel/SEO2/INTEGR~1/variable.genes/ls_genes.Rds")


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
# wbarplot(nrGenesPerDataset, incrBottMarginBy = 2, tilted_text = T
#          , hline = 20000, w = 14, h=7, )

# End ------------------------
create_set_Original_OutDir()



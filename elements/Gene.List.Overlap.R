######################################################################
# Gene.List.Overlap.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/Gene.List.Overlap.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try(source('~/GitHub/Packages/Seurat.utils/Jaccard.toolkit.R'))
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

# idxEx <- which(names(ls_genes) %in% "TSC.122580_WT.20201015")
# if (l(which(names(ls_genes) %in% "TSC.122580_WT.20201015"))) {
#   ls_genes <- ls_genes[-idxEx]
# }


# str(ls_genes)




# Calculate Jaccard ------------------------
df.presence <- jPresenceMatrix(string_list = ls_genes)
PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)

# ------------------------
colnames(PairwiseJaccardIndices)

ph.Jac.wo <- pheatmap::pheatmap(PairwiseJaccardIndices, display_numbers = T, cutree_rows = 3, cutree_cols = 3)
wplot_save_pheatmap(ph.Jac.wo, width = 15)

# End ------------------------
create_set_Original_OutDir()



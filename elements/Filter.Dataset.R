######################################################################
# Filter.Dataset.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/Filter.Dataset.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Filter.Dataset.R")
# try(dev.off(), silent = T)

# Setup ------------------------
# Calculate ------------------------


tic();
# ls.Seurat <- foreach(i = 1:n.datasets ) %dopar% {
for (i in 1:n.datasets ) {
  iprint(names(ls.Seurat)[i], percentage_formatter(i/n.datasets, digitz = 2))
  sobj = ls.Seurat[[i]]
  sobj = subset(x = sobj, subset = `nFeature_RNA` > p$'thr.hp.nFeature_RNA' & `nFeature_RNA` < p$'thr.lp.nFeature_RNA')
  sobj = subset(x = sobj, subset = `percent.mito` > p$'thr.hp.mito' & `percent.mito` < p$'thr.lp.mito')
  sobj = subset(x = sobj, subset = `percent.ribo` > p$'thr.hp.ribo' & `percent.ribo` < p$'thr.lp.ribo')
  ls.Seurat[[i]] <- sobj
}; toc();


Nr.Cells.After.Filtering = unlapply(ls.Seurat, ncol); names(Nr.Cells.After.Filtering) = samples.short

CellCountsAndFiltering <- cbind(
  "Before" = Nr.Cells.Before.Filtering
  , "After" = Nr.Cells.After.Filtering
)

if (n.datasets>1) {
  wplot(CellCountsAndFiltering, abline = "ab", a=0, b=1
        , col = wcolorize(meta.tags$library, set = "rich")
        , xlim = c(0,max(Nr.Cells.Before.Filtering))
        , ylim = c(0, max(Nr.Cells.Before.Filtering))
  )
}

rm('sobj')

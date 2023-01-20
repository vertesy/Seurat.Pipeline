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
  {
    "Dynamic nFeature LP cutoff at 99.75% percentile"
    below.nFeature_RNA <- floor(quantile(ls.Seurat[[i]]$'nFeature_RNA', probs = 0.9975))
  }
  sobj = subset(x = sobj, subset = `nFeature_RNA` > p$'thr.hp.nFeature_RNA' & `nFeature_RNA` < below.nFeature_RNA) # p$'thr.lp.nFeature_RNA'
  sobj = subset(x = sobj, subset = `percent.mito` > p$'thr.hp.mito' & `percent.mito` < p$'thr.lp.mito')
  sobj = subset(x = sobj, subset = `percent.ribo` > p$'thr.hp.ribo' & `percent.ribo` < p$'thr.lp.ribo')
  ls.Seurat[[i]] <- sobj
}; toc();


Nr.Cells.After.Filtering = unlapply(ls.Seurat, ncol); names(Nr.Cells.After.Filtering) = samples.short

CellCountsAndFiltering <- cbind(
  "Before" = Nr.Cells.Before.Filtering
  , "After" = Nr.Cells.After.Filtering
)
rownames(CellCountsAndFiltering) <-  samples.short



if (n.datasets>1) {

  dotlabels <- p0(rownames(CellCountsAndFiltering), '\n'
                  , p0(CellCountsAndFiltering[,'After'], ' cells'), '\n'
                  , percentage_formatter(CellCountsAndFiltering[,'After']/CellCountsAndFiltering[,'Before']
                                         )
    )

  qscatter(CellCountsAndFiltering
           , subtitle = "Fraction of cells remaining after filtering"
           , abline = c(0,1)
           , label = dotlabels
           , xlim = c(0,max(Nr.Cells.Before.Filtering))
           , ylim = c(0, max(Nr.Cells.Before.Filtering)))

  # wplot(CellCountsAndFiltering, abline = "ab", a=0, b=1
  #       , col = wcolorize(meta.tags$library, set = "rich")
  #       , xlim = c(0,max(Nr.Cells.Before.Filtering))
  #       , ylim = c(0, max(Nr.Cells.Before.Filtering))
  # )
}

rm('sobj')


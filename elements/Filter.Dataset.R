######################################################################
# Filter.Dataset.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/Filter.Dataset.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Filter.Dataset.R")
# try(dev.off(), silent = T)

# Setup ------------------------
ScatPlotSize <- max(6, round(1.2 *n.datasets))

# Plot ------------------------

Nr.Cells.Before.Filtering <- unlapply(ls.Seurat, ncol); # names(Nr.Cells.Before.Filtering) = suffices
qbarplot(Nr.Cells.Before.Filtering, label = Nr.Cells.Before.Filtering, ylab = "Cells"
         , subtitle = paste('Total:', sum(Nr.Cells.Before.Filtering)) )


# Filter ------------------------
tic();
for (i in 1:n.datasets ) {
  iprint(names(ls.Seurat)[i], percentage_formatter(i/n.datasets, digitz = 2))
  sobj = ls.Seurat[[i]]

  {
    "Dynamic nFeature LP cutoff at 99.75% percentile"
    stopifnot( !is_null(p$'quantile.thr.lp.nFeature_RNA'))
    stopifnot(p$'quantile.thr.lp.nFeature_RNA'>0 & p$'quantile.thr.lp.nFeature_RNA' <1)
    below.nFeature_RNA <- floor(quantile(ls.Seurat[[i]]$'nFeature_RNA', probs = p$'quantile.thr.lp.nFeature_RNA'))
  }
  sobj = subset(x = sobj, subset = `nFeature_RNA` > p$'thr.hp.nFeature_RNA' & `nFeature_RNA` < below.nFeature_RNA) # p$'thr.lp.nFeature_RNA'
  sobj = subset(x = sobj, subset = `percent.mito` > p$'thr.hp.mito' & `percent.mito` < p$'thr.lp.mito')
  sobj = subset(x = sobj, subset = `percent.ribo` > p$'thr.hp.ribo' & `percent.ribo` < p$'thr.lp.ribo')
  ls.Seurat[[i]] <- sobj
}; toc();


Nr.Cells.After.Filtering <- unlapply(ls.Seurat, ncol)
names(Nr.Cells.After.Filtering) <- samples.short


qbarplot(Nr.Cells.After.Filtering, label = Nr.Cells.After.Filtering, ylab = "Cells"
         , subtitle = paste('Total:', (sum(Nr.Cells.After.Filtering))))

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

  plotmin = round(0.9 * min(Nr.Cells.Before.Filtering))
  plotmax = 1.2 * max(Nr.Cells.Before.Filtering)
  ggExpress::qscatter(CellCountsAndFiltering
           , subtitle = "Fraction of cells remaining after filtering"
           , abline = c(0,1)
           , label = dotlabels
           , xlim = c(plotmin, plotmax)
           , ylim = c(plotmin, plotmax)
           , w = ScatPlotSize, h = ScatPlotSize
           # , col = 1:nrow(CellCountsAndFiltering)
  )


}

rm('sobj')


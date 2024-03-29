######################################################################
# Dowsample.Seurat.Objects.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/Dowsample.Seurat.Objects.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Dowsample.Seurat.Objects.R")


# Metadata ------------------------
# Parameters ------------------------
ls.Seurat.downsampled <- list.fromNames(names(ls.Seurat))
# ls.Seurat.downsampled <- foreach(i=1:n.datasets ) %dopar% {
  for (i in 1:n.datasets ) {
  iprint(names(ls.Seurat)[i], percentage_formatter(i/n.datasets, digitz = 2))
  sobj <- downsampleSeuObj(ls.Seurat[[i]], nCells = p$"dSample.Organoids")
  ls.Seurat.downsampled[[i]] <- sobj
  # sobj
}; toc();# names(ls.Seurat)  <- samples
unlapply(ls.Seurat.downsampled, ncol)
isave.RDS(object = ls.Seurat.downsampled, suffix = ppp(p$"dSample.Organoids", "cells"), inOutDir = T)


# ------------------------
unlapply(ls.Seurat, ncol)
# ------------------------
# ------------------------
# ------------------------
# ------------------------



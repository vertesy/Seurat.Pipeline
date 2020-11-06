######################################################################
# Renumber.Clusters.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/Renumber.Clusters.R")


# Functions ------------------------
source('~/GitHub/Packages/Seurat.utils/Cluster.Auto-naming.DE.R')
stopifnot(exists('p'))
stopif(is.null(p$'OrderClsByPrCurve'),message = 'OrderClsByPrCurve not found')
stopif(is.null(p$'res.analyzed.DE'),message = 'res.analyzed.DE not found')

if (p$'OrderClsByPrCurve') require(princurve)

# Setup ------------------------
create_set_Original_OutDir()
suffix = ""
prefix = p0(if (p$'integrate.multiple') "integrated" else "RNA", "_snn_res")


# AutoNumber.by.UMAP ------------------------
for (i in 1:length(p$'res.analyzed.DE')) {
  res = p$'res.analyzed.DE'[i]
  (res.full = trimws(ppp(prefix, res, suffix), whitespace = '\\.'))

  combined.obj <- AutoNumber.by.UMAP(obj = combined.obj, dimension=1, reduction="umap", swap = F, res = res.full )
  identity.used <- ppp(res.full, "ordered")
  umap.snn_res.X.ordered <- DimPlot.ClusterNames(ident = identity.used)
  qqsave(umap.snn_res.X.ordered
         , plotname = ppp("umap",identity.used), PNG = T, plotit=T)


  # AutoNumber.by.PrinCurve ------------------------
  if (p$'OrderClsByPrCurve') {
    combined.obj <- AutoNumber.by.PrinCurve(obj = combined.obj, dimension=1:2, reduction="umap", plotit=T, swap= -1, res = res.full )
    identity.used <- ppp(res.full, "prin.curve")
    umap.snn_res.X.prin.curve <- DimPlot.ClusterNames(ident = identity.used)
    qqsave(umap.snn_res.X.prin.curve
           , plotname = ppp("umap",identity.used), PNG = T, plotit=T)
  }


}




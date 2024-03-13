######################################################################
# Renumber.Clusters.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/Renumber.Clusters.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Renumber.Clusters.R")


# Functions ------------------------
# source('~/GitHub/Packages/Seurat.utils/Functions/Cluster.Auto-naming.DE.R')
stopifnot(exists('p'))
stopif(is.null(p$'OrderClsByPrCurve'), message = 'OrderClsByPrCurve not found')
stopif(is.null(p$'res.analyzed.DE'), message = 'res.analyzed.DE not found')

if (p$'OrderClsByPrCurve') require(princurve)

# Setup ------------------------
create_set_Original_OutDir()
suffix = ""
if(!exists('n.datasets')) n.datasets <- 1
(prefix = p0(if (n.datasets > 1 & p$'integration' == 'CCA') "integrated" else "RNA", "_snn_res"))


# AutoNumber.by.UMAP ------------------------
for (i in 1:length(p$'res.analyzed.DE')) {
  res = p$'res.analyzed.DE'[i]
  (res.full = trimws(ppp(prefix, res, suffix), whitespace = '\\.'))

  combined.obj <- AutoNumber.by.UMAP(obj = combined.obj, dim =  abs(p$"Reorder.Dim"), reduction = "umap", swap = (p$"Reorder.Dim" < 0), res = res.full)
  identity.used <- ppp(res.full, "ordered")
  clUMAP(identity.used)

  # AutoNumber.by.PrinCurve ------------------------
  if (p$'OrderClsByPrCurve') {
    combined.obj <- AutoNumber.by.PrinCurve(obj = combined.obj, dim = 1:2, reduction = "umap", plotit = T, swap = (p$"Reorder.Dim" < 0), res = res.full)
    identity.used <- ppp(res.full, "prin.curve")
    clUMAP(identity.used)
  }

}




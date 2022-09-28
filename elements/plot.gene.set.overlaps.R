######################################################################
# plot.gene.set.overlaps.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/plot.gene.set.overlaps.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Plot.gene.set.overlaps.R")
# Based on https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
# try(dev.off(), silent = T)

# Parameters ------------------------
p$"plotUpset" <- F
# p$"n.overlaps" = 5

# Plot ------------------------
ls_genes = lapply(ls.Seurat, rownames)
# names(ls_genes) = samples.short
try.dev.off()
ls_genes = lapply(ls.Seurat, rownames);
names(ls_genes) = samples.short; try.dev.off()

gene.overlap <- UpSetR::fromList(ls_genes)

if (length(ls_genes) < 6) {
  # Wenn diagram ------------------------------------------------
  wvenn(ls_genes)

  # plotUpset ------------------------------------------------
  # devtools::install_github("hms-dbmi/UpSetR") # too old version: install.packages("UpSetR")

} else if (p$"plotUpset" & (length(ls_genes) < 12) ) {
  library("UpSetR")

  # Sys.setenv('R_MAX_VSIZE'=32000000000)
  p$"n.overlaps" = 5
  upsetplot <- upset(data = gene.overlap,  sets.bar.color = "#56B4E9"
                     # , boxplot.summary = T
                     , nsets = l(ls_genes), nintersects = p$"n.overlaps"
                     , empty.intersections = "on" , order.by = "freq"
  )
  upsetplot
  # wplot_save_this(plotname = 'upsetplot')
  ggsave2(filename = "Overlaps.png" , plot = upsetplot)

} else {
  iprint("You need to run Gene.List.Overlap.R")
}

# barplot ------------------------------------------------
rowSums(gene.overlap)
GenesDetected <- colSums(gene.overlap)
qbarplot(GenesDetected, xlab.angle = 45)
# wbarplot(GenesDetected, incrBottMarginBy = 3, tilted_text = 15
#          , col = wcolorize(meta.tags$project, set = "rich"))


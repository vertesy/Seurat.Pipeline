######################################################################
# STRING.db.auto.cluster.annotation.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/STRING.db.auto.cluster.annotation.R')
library(magrittr)

# Setup ------------------------
create_set_OutDir(OutDirOrig, "STRING.db.auto.cluster.annotation")

if (F) {
  # install.packages("BiocManager")
  if (!requireNamespace("BiocManager", quietly = TRUE))
  BiocManager::install("STRINGdb")
} else { irequire("STRINGdb") }


# Parameters ------------------------
GO.terms.removed = c("extracellular structure organization")

GetNamedClusteringRuns
p$'Ident.for.STRING' <- p$'Ident.for.DEG'
# gsub(p$'Ident.for.DEG', pattern = '[0-9]\\.[0-9]'
#                               , replacement = 'STRING.res', perl = T)

GOident = ppp(GetClusteringRuns(res = p$def_res), 'GO.process')

# STRING database ----------------------
Idents(combined.obj) <- p$'Ident.for.DEG'
clUMAP(ident = p$'Ident.for.DEG')

v.clusters <- unique(Idents(combined.obj))
string_link <- vec.fromNames(v.clusters)
DEG.top.genes.Padj <- list.fromNames(v.clusters)

# x=get_STRING_species(version = 11, species_name=NULL); View(x)
# Homo sapiens 9606
# Mus musculus 10090
if (!exists("string_db")) string_db <- STRINGdb$new(version = "11", species = 9606, score_threshold = 0, input_directory = "" )


i = 15

clnames.GO <- vec.fromNames(sort(v.clusters))
clzUsed <- (as.numeric(as.character(sort(v.clusters))))
if (min(clzUsed) < 1) clzUsed <- clzUsed + 1
for (i in clzUsed) { print(i)
  cl <- (i - 1)
  cl.char <- as.character(cl)

  # table(DEG.top.genes.Padj[[i]]$cluster)
  DEG.top.genes.Padj[[i]] <-
    combined.obj@misc$df.markers[[p$'def_res']] %>%
    # rownames_to_column('gene') %>%
    arrange(p_val_adj) %>%
    filter(cluster == cl) %>%
    filter(avg_logFC > 0.5) %>%     # which means at least 2 fold enriched (e logN()s)
    filter(p_val_adj < 0.05) %>%
    head(n = p$'STRING.nr.genes') %>%
    pull("gene" )

  (genesX = DEG.top.genes.Padj[[i]])
  # write_clip(genesX)
  iprint("    ",l(genesX),"Genes met enrichment criteria.")

  string_ids = string_db$mp(genesX)
  if (l(string_ids) == 0) {
    print("string_ids  is empty")
    next
  }

  if (p$"getSTRINGlinks") {
    string_link[i] = string_db$get_link(string_ids)
    browseURL(string_link[i])
  }

  pl <- qUMAP(genesX[1]) +  ggtitle(label = genesX[1], subtitle = paste('Top marker for cluster',cl))
  save_plot(plot = pl, filename = ppp('STRING.top.gene.cl',cl,genesX[1],'png'))

  # Automatic cluster annotation ------------------------------------------------
  # http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_BS32010/Workshops/BS32010Workshop4.html
  enrichmentGO <- string_db$get_enrichment(string_ids, category = "Process") # , methodMT =  "FDR", iea = TRUE
  enrichmentGO <- enrichmentGO[order(  enrichmentGO$'p_value'),]
  specific.GO.terms <- which(enrichmentGO$number_of_genes_in_background < p$'STRING.max.genes')

  filter.specific.GO.terms = T #
  if (filter.specific.GO.terms) {
    idx.exclude <- which(enrichmentGO$description %in% GO.terms.removed)
    if (l(idx.exclude) > 0) enrichmentGO <- enrichmentGO[-idx.exclude,]
  }

  clnames.GO[cl.char] <-
    if (nrow(enrichmentGO)) {
      if (l(specific.GO.terms)) {  enrichmentGO$description[specific.GO.terms][1]
      } else {                     enrichmentGO$description[1] }
    } else {ppp('No Enrichment', i)} # if

  (clnames.GO[cl.char] <- ppd(clnames.GO[cl.char],l(DEG.top.genes.Padj[[i]])))
}
# write_clip(enrichmentGO$description[specific.GO.terms])


if (p$"getSTRINGlinks") write.simple.tsv(sort.natural(string_link))
clnames.GO <- ppp(v.clusters,clnames.GO)
is(Idents(combined.obj))
symdiff(Idents(combined.obj),  v.clusters)


(Ident.GO <- translate(vec = as.character(combined.obj@meta.data[, p$'Ident.for.DEG'])
                       , oldvalues = as.character(v.clusters), newvalues = clnames.GO))
combined.obj[[GOident]] <- Ident.GO
clUMAP(ident = GOident , sub = "nr. of genes provided to STRING after dash")
say()

# # Manual annotation ------------------------------------------------
# if (F) {
#   clusterIDs.GO.process <- read.simple.tsv('~/Data/POL/metadata/Manual.Annotation.GO.Process.2020.08.26.txt')
#   clusterIDs.GO.process$`GO-Process` <- gsub("\\s*\\([^\\)]+\\)","",clusterIDs.GO.process$`GO-Process`)
#
#   combined.obj <- seu.map.and.add.new.ident.to.meta(obj = combined.obj, ident.table = clusterIDs.GO.process
#                                                     , metaD.colname = 'GO.TOP.process.r.0.5')
#
#   unique(combined.obj@meta.data$GO.TOP.process.r.0.5)
#   clUMAP(ident = 'GO.TOP.process.r.0.5')
# }
# combined.obj@meta.data$integrated_snn_res.0.5.ordered
#
# # clUMAP(ident = 'integrated_snn_res.0.5.ordered', h = 5)
# qUMAP("log10.HGA_Markers")
# qUMAP("S100B")



# Manual annotation ------------------------------------------------
# "select.cells not defined"
# if (condition) {
#   for (ds in (samples.short)) {
#     cells.in.DS <- getCellIDs.from.meta(obj = combined.obj, ColName.meta = 'sample', values = ds)
#
#
#     # Write out for MSeq ----
#     write.simple.vec(input_vec = substr(cells.in.DS, 1, 16)
#                      , ManualName = kpp('CBC.Seurat.filtered', ds,'tsv'))
#
#     # Pot QC ----
#     ccc <- 1+as.numeric(cells.in.DS %in% select.cells)
#     table(ccc)
#     p <- FeatureScatter(object = combined.obj, feature1 = "percent.ribo", feature2 = "percent.mito"
#                         , cells =cells.in.DS
#                         # , cols = ccc
#                         , slot = 'counts') +
#       scale_x_log10() + scale_y_log10() + annotation_logticks(sides="trbl") + ggtitle("Quality vs. Clusters")
#     ggsave(p, filename = kpp("Quality vs. Clusters",ds,"pdf"), width = 7, height = 5)
#   }
# }


# End ------------------------------------------------------------------------
create_set_Original_OutDir()



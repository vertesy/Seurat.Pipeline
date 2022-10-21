######################################################################
# Visualize single-cell gene expression
######################################################################
# source ('')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)
stop()

# Custom functions you need--------------------------------------------------------------------------------
source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/Load.packages.local.R')

# Setup ------------------------
"Your files will be saved under the documents folder"

OutDirOrig = OutDir = p0("~/single-cell-visualization/Analysis/", ppp("Results",idate()))
iprint("Your results will be saved in:", OutDirOrig)

setup_MarkdownReports(OutDir = OutDir, scriptname = "Visualize.Gene.expression.R")

# Read In the experiment you want ------------------------
combined.obj = readr::read_rds('Path/to/the/Seurat_object.Rds.gz')

"See whats it:"
combined.obj

# Show clusters --------------------------------------------------------------------------------
"Visualize  UMAP (better for dev. trajectories)"


(clusterID = GetNamedClusteringRuns(obj = combined.obj)[1])
# You can also try:
clusterID2 = GetClusteringRuns(obj = combined.obj)[1]

clUMAP(ident = clusterID, obj = combined.obj)

"for help on function arguments use ? before comandname"
"See:"
?ggplot
?clUMAP




# Show gene expression --------------------------------------------------------------------------------

all.genes = list.fromNames(rownames(combined.obj), fill = 0)
"find genes by all.genes$START_TYPING and hit tabulator"
"all.genes$start_typing..."


qUMAP(feature = "TOP2A", obj = combined.obj)

"See saved plot in 'OutDir': "
oo()


# Show key marker genes --------------------------------------------------------------------------------
qMarkerCheck.BrainOrg(obj = combined.obj)


# Plot custom gene sets  --------------------------------------------------------------------------------

SEGA.old = c("HES5", "ID3", "STMN1", "HMGN2", "SLC6A9", "VCAM1", "S1PR1", "ETNPPL")

multiFeaturePlot.A4(list.of.genes = SEGA.old, obj = combined.obj)



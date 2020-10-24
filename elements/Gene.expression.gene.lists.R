######################################################################
# Gene.expression.gene.lists.R
######################################################################
# source ('~/GitHub/Packages/Seurat.pipeline/elements/Gene.expression.gene.lists.R')


# create_set_SubDir("Gene.expression")
create_set_OutDir(OutDirOrig,"Gene.expression")
# OutDirOrig ="~/Dropbox/Abel.IMBA/AnalysisD/Bajaj/iN.migration..premRNA."


# Plot QC gene sets -----------------------------------

Highest.Expressed.Genes = names(head(sort(combined.obj@misc$expr.q90, decreasing = T), n = 16))
multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, obj = combined.obj, intersectionAssay = "RNA", subdir =T)

AutoNaming.Genes = p$Cluster.Labels.Hybrid.Genes
multiFeaturePlot.A4(list.of.genes = AutoNaming.Genes, obj = combined.obj, intersectionAssay = "RNA", subdir =T)


# Classic markers ------------------------------------------------------------------------
ClassicMarkers = c(
  "Apical precursor (Dorsal)" =		"SOX2",
  "Stem cells" = 						  		"ID4",
  "Cycling cells" = 							"TOP2A",
  "Intermediate progenitor" = 		"EOMES",
  "Intermediate progenitor1" = 		"TAC3",
  "Immature neurons" = 		 		 		"SLA",

  "adult RGL/NSC at dSVZ, astr.lin" ="HOPX", # https://doi.org/10.1016/j.stemcr.2018.08.006

  "Late Progenitor" = 						"EGFR",
  "Astroglia" = 									"GFAP",
  "Astrocyte" = 									"S100B",
  "Interneurons, (OPC)"   = 			"DLX6-AS1",
  "Interneurons"   = 							"DLX5",
  "Interneurons"   = 							"DLX2",
  "Interneurons, (SST)"   = 			"SST",
  "GABA-ergic inter neuron" = 		"GAD2",

  'Hypoxia/Stress'=   	   	   	  "DDIT4",
  "Glycolytic" =   								"PDK1",

  "CTIP2 " = 	  									"BCL11B",
  "SATB2 " =   										"SATB2",
  "CTIP2" = 											"FEZF2"
)

QC.markers = c(
  "MALAT1" = 										"MALAT1",
  "RPL34" = 										"RPL34",
  "Mito" =   										"MT-ATP6",
  "Glycolysis" =  							"PDK1",
  "IGFBP5" = 										"IGFBP5"
)

LargeSubsetMarkers = c(

  "Dorsal neuronal" = 			"NEUROD6",
  "non-neuronal cells" = 		"VIM", # https://www.labome.com/method/Neuronal-Cell-Markers.html
  "All neurons" = 	      	"DCX",
  "All neurons" =          	'MAP2',
  "All progenitors" = 			'SOX6',


  'Dorsal neurons' = 				'RBFOX3',
  "Dorsal progenitors" = 			'PAX6',
  "Dorsal progenitors" = 			'SOX3',
  "Dorsal progenitors" = 			'FOS',
  'Dorsal progenitors' = 		'FOXP2')


# ClassicMarkers.found <- check.genes(obj = combined.obj, list.of.genes = ClassicMarkers)
ClassicMarkers.plus = c(ClassicMarkers, LargeSubsetMarkers, QC.markers)
multiFeaturePlot.A4(list.of.genes = ClassicMarkers.plus, obj = combined.obj, subdir =T)
# multiFeaturePlot.A4(list.of.genes = ClassicMarkers, obj = combined.obj, plot.reduction = 'tsne', subdir =T)

AdditionalMarkers = c("BRN2", "TBR1", "TBR2", "VGLUT1", "VGLUT2", "NCAM", "L1CAM", "PSD95")
multiFeaturePlot.A4(list.of.genes = AdditionalMarkers, obj = combined.obj, subdir =T)


# End ------------------------------------------------------------------------
create_set_Original_OutDir()


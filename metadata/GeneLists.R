######################################################################
# GeneLists.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/metadata/GeneLists.R')


# Metadata ---------------------------
HGA_MarkerGenes <- c("ENO1", "IGFBP2", "WSB1", "DDIT4", "PGK1", "BNIP3", "FAM162A",
                     "TPI1", "VEGFA", "PDK1", "PGAM1", "IER2", "FOS", "BTG1", "EPB41L4A-AS1",
                     "NPAS4", "HK2", "BNIP3L", "JUN", "ENO2", "GAPDH", "ANKRD37",
                     "ALDOA", "GADD45G", "TXNIP")

# Classic markers ------------------------------------------------------------------------
ClassicMarkers = c(
  "Apical precursor (Dorsal)" =		"SOX2",
  "Stem cells" = 						  		"ID4",
  "Cycling cells" = 							"TOP2A",
  "Intermediate progenitor" = 		"EOMES",
  "Intermediate progenitor1" = 		"TAC3",
  "Immature neurons" = 		 		 		"SLA",
  
  "adult RGL/NSC at dSVZ, astr.lin" = "HOPX", # https://doi.org/10.1016/j.stemcr.2018.08.006
  
  "Late Progenitor" = 						"EGFR",
  "Astroglia" = 									"GFAP",
  "Astrocyte" = 									"S100B",
  "Interneurons, (OPC)"   = 			"DLX6-AS1",
  "Interneurons"   = 							"DLX5",
  "Interneurons"   = 							"DLX2",
  "Interneurons, (SST)"   = 			"SST",
  "GABA-ergic inter neuron" = 		"GAD2",
  
  'Hypoxia/Stress' =   	   	   	  "DDIT4",
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
  'Dorsal progenitors' = 		'FOXP2'
)



######################################################################
# Parameters.example.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/Parameters.example.R')

'Not all parameters are used'
# Parameters ------------------------
p = NULL

{
  # Data and Script parameters (subscript settings) -----
  p$"file.ext"                      = ".png"
  p$'parallelize.w.futures'         = F
  p$'premRNA'                       = T           # Did you use premRNA transcript model in mapping? Used for folder naming.
  p$'MemSaving'                     = T           # Remove unused data objects from memory.
  p$'seed'                          = 1989
  p$"Filter.HGA.clusters"           = F           # create cell whitelists based on exclusion of glycolytic cells
  p$'Subset'                        = T           # Subsetting or pre-Filtering for non-stressed cells.
  p$'cellWhiteList'                 = '~/Dropbox/Abel.IMBA/MetadataD/INM.meta/cell.lists/' # List of cells taken for the analysis in a 2 column table: (cellID - sampleID)
  p$'subclustering'                 = F           # Is this run sub-clustering run?

  # Filtering parameters -----
  # LP for low-pass (below), HP for high-pass (above)
  p$'min.cells'             = 10
  p$'min.features'          = 200
  p$'thr.lp.nFeature_RNA'   = 6000
  p$'thr.hp.nFeature_RNA'   = 900
  p$'thr.hp.mito'           = 0.01                  # Minimal fractions of mitochondrial reads in total transcriptome.
  p$'thr.lp.mito'           = 0.15
  p$'thr.hp.ribo'           = 0.05
  p$'thr.lp.ribo'           = 0.4

  p$"res.stress.filtering"  = 0.5

  # Seurat parameters -----
  p$'dSample.Organoids'     = F                     # Do you want to downsample each the dataset (unmerged Seurat object) to a given number of cells?
  p$'variables.2.regress'   = NULL                  # c( 'percent.mito', 'percent.ribo', 'nFeature_RNA')
  p$'n.CC'                  = 50                    # How many CCA dimensions to use to specify the neighbor search space [FindIntegrationAnchors]
  p$'n.PC'                  = p$'n.CC'              # How many PCA dimensions to use for UMAP & co
  p$"flipUmapCoord"         = 'x'                   # Do you want to Flip x/y UMAP coordinates?
  p$'snn_res'               = c(.3, .5, .7, .8)     # Which clustering resolutions to calculate
  p$'def_res'               = .7                    # Which clustering resolutions to calculate

  p$'res.analyzed.DE'       = c(.3, p$'def_res')    # Which clustering resolutions to analyze in the differential gene expression script
  p$'n.markers'             = 7                     # The top how many markers to plot per cluster
  # p$'integration'          = 'CCA'                # What integration method to use?



  # Plot parameters - basic plots -----
  p$"plotHexBinStatPlots"     = T
  p$"plotClusterPhylogeny"    = T
  p$"StatFeatures"            = c( "percent.mito", "percent.ribo", "nFeature_RNA", "nCount_RNA")

  # Differential Gene Expression parameters -----
  # See ?FindAllMarkers

  p$"cl.annotation"    	  	= c("ordered", "named", "simple")[1] #
  p$"only.pos" 		    	  	= TRUE
  p$"test"     		    	  	= c("wilcox", "MAST")[1] # Two best tests
  p$"min.pct" 				    	= 0.1               # Min 10% cells need to be expressing it.
  p$"min.diff.pct" 					= 0.01              # Def 0 | Removes <100 very prevalent genes from DE results. Some, like GNAS, is clearly DE, but its not a unique marker for sure
  p$"min.cells.group" 			= 20                # Def 3 | Minimum number of cells in one of the groups
  p$"min.cells.feature" 		= 20                # Def 3 | Minimum number of cells expressing the feature in at least one of the two groups, currently only used for poisson and negative binomial tests
  p$"logfc.threshold" 			= 0.25              # Def 0.25 | Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.
  p$"return.thresh"     		= 1e-3              # Def 1e-2 | Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
  p$"DEG.ranking"        		= c("combined.score","avg_logFC", "p_val_adj")[1]      #


  # Cluster Annotation -----
  p$'Cluster.Labels.Automatic'      = T    # Label clusters by top DE gene

  p$'Cluster.Labels.Hybrid.Genes'   = c(" G2M-phase" = 'TOP2A' # , 'MKI67'
                                        , "S-phase" = 'HIST1H4C' # S phase
                                        , 'NES' # neuron fate commitment. Also DLL3+ cells
                                        , "Migrating" = 'CXCR4' # migrating https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025760/
                                        , 'SCGN' # CXCR4+ cells too
                                        , 'GAD1'
                                        , 'MEF2C' # Gene list to label clusters by top DE gene intersected with a predefined set of genes
                                        , "Serotonin" = 'MEIS2'
                                        , "DLX6"
                                        , "LHX8"
                                        , 'NR2F2'
                                        # , 'ZFHX3'
                                        , "Somatostatin iN" = 'SST'
                                        , "Pallial iN" = 'MAF'
                                        # , 'ERBB4' # ALSO MAF+
                                        # , 'SIX3'  # Medium Spiny Neuron https://dev.biologists.org/content/145/14/dev165456.long OR cholinergic interneurons https://europepmc.org/article/pmc/pmc2786914
                                        , 'NKX2-1'
                                        # , 'PLS3' # Subplate https://www.sciencedirect.com/science/article/pii/S1934590917302862
                                        # , 'EBF1' # Striatal neurons that are also 'FOXP1' + https://dev.biologists.org/content/126/23/5285.long
                                        # , 'DNAJA1' # Heatshock cluster
                                        , 'NFIA' # Glia https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6591152/
                                        , 'FOXP1' # Projection
                                        # , 'ISL1' # Also FOXP1
                                        )

  # qUMAP('LHX8')
  # clUMAP("cl.names.KnownMarkers.0.7")

  # 3D plotting -----
  p$"cex3D"                     	= 3
  p$"CustomGenesFor3Dplots"      	= c("CALB1", "CALB2", "CCND2", "CDK14", "DLX1", "DLX2", "DLX5",
                                "DLX6", "EBF1", "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5",
                                "EFNB1", "EFNB2", "EFNB3", "EOMES", "EPHA4", "EPHA5", "ERBB4",
                                "FOXP1", "FOXP2", "GAD1", "GAD2", "GFAP", "GSX2", "HTR3A", "ID2",
                                "IGFBP6", "ISL1", "LHX6", "MAF", "MAP2", "MEF2C", "MEIS2", "MKI67",
                                "MYT1L", "NEUROD6", "NKX2-1", "NOS1", "NPY", "NR2F2", "OLIG1",
                                "OLIG2", "PDGFRA", "PLCXD3", "PROX1", "PVALB", "RELN", "S100B",
                                "SATB1", "SCGN", "SLC17A6", "SLC17A7", "SLC32A1", "SOX2", "SOX6",
                                "SP8", "SST", "ST18", "STMN2", "TBR1", "TH", "TSPAN7", "VIM",
                                "VIP", "ZFHX3")


  # p$'Cluster.Labels.Hybrid.Genes'   = c('CALB2', 'ERBB4', 'FOXP1', 'GAD1', 'ID2', 'ISL1', 'MAF', 'MEF2C', # Gene list to label clusters by top DE gene intersected with a predefined set of genes
  #                                       'MEIS2', 'MKI67', 'NKX2-1', 'NR2F2', 'SCGN', 'SST', 'ZFHX3')
  # For normal DF organoids: c('TOP2A', 'EOMES', 'SLA', 'HOPX', 'S100B', 'DLX6-AS1', 'POU5F1','SALL4','DDIT4', 'PDK1', 'SATB2', 'FEZF2')
  p$'auto.name.res'                 = p$'def_res'             # Which clustering resolution to use for to autoname clusters (e.g. by top gene)?
  p$'auto.name.method'              = c('topgene', 'hybrid.intersect', 'hybrid.identify')[1:2]    # What cluster autonaming methods to perform? E.g.: can be named by the top gene.
  p$'res.MetaD.colname'             = ppp("cl.names.KnownMarkers",p$'def_res')
  # p$'res.MetaD.colname'             = GetNamedClusteringRuns(res = p$'auto.name.res')
  p$"Reorder.Dim"                   = c(1, -1, 2, -2)[1] # By which UMAP dimension, and in which direction should renumbering take place?
  p$'OrderClsByPrCurve'             = F    # Re-order clusters numbers by fitting a principal curve

  p$'STRING.nr.genes'             = 100
  p$"getSTRINGlinks"                = F
  assay.def = 'RNA'

  # Other less used parameters -----
  # p$'Man.Int.Order'           = F
  # p$'with.Nowakowski'         = F
  # p$'Update'                  = F                     # Do you want to update the Seurat object?
  # p$'use.SCTransform'               = F # Not implemented - worked very badly
  # p$'with.ORC'                      = F    # Analyze with Bhaduri
  # p$'Cluster.Labels.Manual'         = T
  # p$'auto.name.by.classic.markers'  = T               # Not used. Control if auto.name.by.classic.markers is done

}

# Checks ----------------------------------------------------------------------------------------------------
stopifnot(
  all(
    p$'auto.name.res' %in% p$'snn_res',
    all(p$'res.analyzed.DE' %in% p$'snn_res')
  )
)


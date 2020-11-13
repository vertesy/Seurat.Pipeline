######################################################################
# JackStraw.plot.R
######################################################################
# source('~/GitHub/Packages/Seurat.pipeline/elements/JackStraw.plot.R') # takes very long
try.dev.off()
if (F) source('~/GitHub/Packages/Seurat.pipeline/elements/JackStraw.plot.R') # takes very long

# ------------------------
if (F) {
  combined.obj <- readRDS("~/Dropbox/Abel.IMBA/AnalysisD/Abel/CON/scPilot.2.w.RabV.NoFilt.premRNA./00.CON.pilot.2.Frame.2020.11.R_2020_11_13-13h/combined.obj__2020.11.13_13.40.Rds.gz")
  OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/Abel/CON/scPilot.2.w.RabV.NoFilt.premRNA./00.CON.pilot.2.Frame.2020.11.R_2020_11_13-13h/"
}

combined.obj <- JackStraw(combined.obj, reduction = "pca", assay = NULL, dims = 5
                          , num.replicate = 100, prop.freq = 0.01, verbose = TRUE, maxit = 1000)

JackStrawPlot(object = combined.obj, reduction = "pca", dims = 1:5)


######################################################################
# Filtering.plots.3D.solo.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/Filtering.plots.3D.solo.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Filtering.plots.3D.solo.R")
# try(dev.off(), silent = T)



# Setup ------------------------
create_set_OutDir(OutDirOrig, "Filtering.plots3D")
plot3DFiltLinear = F

# stopifnot(length(combined.obj$log10.HGA_Markers) > 0 )


# Calculate ------------------------

for (i in 1:n.datasets) {
  iprint(names(ls.Seurat)[i], percentage_formatter(i/n.datasets, digitz = 1))

  plotting.data <- ls.Seurat[[i]]@meta.data
  plotting.data$'percent.ribo' <- clip.outliers(plotting.data$percent.ribo, probs = c(.005, .995))
  plotting.data$'percent.mito' <- clip.outliers(plotting.data$percent.mito, probs = c(.005, .995))
  plotting.data$'nFeature_RNA' <- clip.outliers(plotting.data$nFeature_RNA, probs = c(.005, .995))

  plotting.data$'percent.ribo.log10' <- log10(plotting.data$percent.ribo + 0.01)
  plotting.data$'percent.mito.log10' <- log10(plotting.data$percent.mito + 0.01)
  plotting.data$'log_nFeature_RNA' <- log10(plotting.data$nFeature_RNA + 1)

  # Plot3D log scale ------------------------
  pltLog <-  plot_ly(data = plotting.data
                 , x = ~log_nFeature_RNA, y = ~percent.mito.log10, z = ~percent.ribo.log10
                 , type = "scatter3d"
                 , mode = "markers"
                 , marker = list(size = 2)
                 , color = ~log10.HGA_Markers
                 , colors = c('blue', 'red')
                 , opacity = .9
                 #, hoverinfo="text"
                 # , text=~label
  ) %>% layout(scene = list(aspectmode = "manual"
                            , title = "Filtering in 3 parameters separate dying cells"
                            , aspectratio = list(x = 1, y = 1, z = 1)))

  # pltLog

  SavePlotlyAsHtml(pltLog, category. = ppp("log10.Filtering", meta.tags$library[i]))

  # Plot3D ------------------------
  if (plot3DFiltLinear) {
    pltLin <-  plot_ly(data = plotting.data
                    , x = ~nFeature_RNA, y = ~percent.mito, z = ~percent.ribo
                    , type = "scatter3d"
                    , mode = "markers"
                    , marker = list(size = 2)
                    , color = ~log10.HGA_Markers
                    , colors = c('blue', 'red')
                    , opacity = .9
                    #, hoverinfo="text"
                    # , text=~label
    ) %>% layout(scene = list(aspectmode = "manual", title = "Filtering in 3 parameters separate dying cells"
                              , aspectratio = list(x = 1, y = 1, z = 1)))
    # pltLin
    SavePlotlyAsHtml(pltLin, category. = ppp("Filtering", meta.tags$library[i]))
  }

}



# End ------------------------------------------------------------------------
create_set_Original_OutDir()

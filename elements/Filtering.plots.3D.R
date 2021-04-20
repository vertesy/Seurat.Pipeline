######################################################################
# Filtering.plots.3D.R
######################################################################
# source("~/GitHub/Packages/Seurat.pipeline/elements/Filtering.plots.3D.R")
# source("https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Filtering.plots.3D.R")
# try(dev.off(), silent = T)

# Setup ------------------------
create_set_OutDir(OutDirOrig, "Filtering.plots.3D")


for (i in 1:n.datasets ) {
  dName <- names(ls.Seurat)[i]
  iprint(dName, percentage_formatter(i/n.datasets, digitz = 2))
  if( !length(ls.Seurat[[i]]$log10.HGA_Markers) > 0 ) { iprint("log10.HGA_Markers is not present"); break  }
  # Calculate ------------------------

  plotting.data <- ls.Seurat[[i]]@meta.data
  plotting.data$'percent.ribo' <- clip.outliers(plotting.data$percent.ribo, probs = c(.005, .995))
  plotting.data$'percent.mito' <- clip.outliers(plotting.data$percent.mito, probs = c(.005, .995))
  plotting.data$'nFeature_RNA' <- clip.outliers(plotting.data$nFeature_RNA, probs = c(.005, .995))

  plotting.data$'percent.ribo.log10' <- log10(plotting.data$percent.ribo + 0.01)
  plotting.data$'percent.mito.log10' <- log10(plotting.data$percent.mito + 0.01)
  plotting.data$'log_nFeature_RNA' <- log10(plotting.data$nFeature_RNA+1)


  # Plot3D ------------------------
  plt <-  plot_ly(data = plotting.data
                  , x = ~log_nFeature_RNA, y = ~percent.mito.log10, z = ~percent.ribo.log10
                  , type = "scatter3d"
                  , mode = "markers"
                  , marker = list(size = 2)
                  , color = ~log10.HGA_Markers
                  , colors = c('blue', 'red')
                  , opacity = .9
                  #, hoverinfo="text"
                  # , text=~label
  ) %>% layout(scene = list(aspectmode = "manual", title="Filtering in 3 parameters separate dying cells"
                            , aspectratio = list(x=1, y=1, z=1)))

  plt
  SavePlotlyAsHtml(plt, category. = ppp(dName, "Filtering.log10"))

  # Plot3Dlog scale ------------------------
  plt <-  plot_ly(data = plotting.data
                  , x = ~nFeature_RNA, y = ~percent.mito, z = ~percent.ribo
                  , type = "scatter3d"
                  , mode = "markers"
                  , marker = list(size = 2)
                  , color = ~log10.HGA_Markers
                  , colors = c('blue', 'red')
                  , opacity = .9
                  #, hoverinfo="text"
                  # , text=~label
  ) %>% layout(scene = list(aspectmode = "manual", title="Filtering in 3 parameters separate dying cells"
                            , aspectratio = list(x=1, y=1, z=1)))
  plt
  SavePlotlyAsHtml(plt, category. = ppp(dName, "Filtering"))


}


# End ------------------------------------------------------------------------
create_set_Original_OutDir()

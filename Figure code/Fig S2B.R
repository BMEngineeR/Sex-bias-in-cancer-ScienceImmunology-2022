rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)

# load data

load( file="seurat_data_aggregated_with_annotation.RData" )

# Visualization

gene.list <- c( "Sell", "Tcf7", "Cd44", "Cd69",
  "Gzma", "Gzmb", "Prf1", "Lag3", 
  "Havcr2", "Pdcd1", "Ctla4", "Tox", 
  "Tbet", "Eomes", "Nr4a1", "Nfatc1",
  "Zeb1", "Slamf6", "Il2ra" )

pdf( "featureplot_markers_full.pdf" )
for ( i in 1:length(gene.list) ) {
  p1 <- FeaturePlot( experiment.aggregate, gene.list[i], col=c("lightyellow","orange","darkred"), min.cutoff="q5", max.cutoff="q90" )
  print(p1)
}
dev.off()

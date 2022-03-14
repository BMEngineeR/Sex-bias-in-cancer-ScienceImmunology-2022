rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)

# load data

load( file="seurat_data_aggregated_with_annotation.RData" )
all.genes <- rownames(experiment.aggregate)
experiment.aggregate <- ScaleData(experiment.aggregate, features = all.genes)

# Fig 3A

print(DimPlot(experiment.aggregate, reduction = "umap", split.by = "sex", label=TRUE ))

# Fig 3B

table( experiment.aggregate@meta.data$sex )
table( experiment.aggregate@meta.data$seurat_clusters )
clust.M <- experiment.aggregate@meta.data$seurat_clusters[ experiment.aggregate@meta.data$sex == "Male" ]
clust.F <- experiment.aggregate@meta.data$seurat_clusters[ experiment.aggregate@meta.data$sex == "Female" ]
( rbind( round( 10000 * table(clust.M) / length(clust.M) ) / 100, 
  round( 10000 * table(clust.F) / length(clust.F) ) / 100 ) )[ ,1:11 ]

# Identify conserved cell type markers

DefaultAssay(experiment.aggregate) <- "RNA"
Idents(experiment.aggregate) <- "seurat_clusters"
list.markers <- mclapply( 1:11, function(i) {
  FindConservedMarkers(experiment.aggregate, ident.1 = (i-1), grouping.var = "sex", verbose = FALSE)  
}, mc.cores=11 )

save( list.markers, file="seurat_list_markers.RData" )
for ( i in 1:length(list.markers) ) {
  write.csv( list.markers[[i]], file=paste("seurat_list_markers_",(i-1),".csv",sep="") )
}

# Identify differential expressed genes between conditions (for Fig 3C and S2A)

experiment.aggregate$celltype.sex <- paste(Idents(experiment.aggregate), experiment.aggregate$sex, sep = "_")
Idents(experiment.aggregate) <- "celltype.sex"
cluster.sex.diff <- vector( "list", 11 )
for ( i in 1:11 ) {
  cluster.sex.diff[[i]] <- FindMarkers(experiment.aggregate, ident.1 = paste((i-1),"_Female",sep=""), ident.2 = paste((i-1),"_Male",sep=""), verbose = FALSE)
}
for ( i in 1:11 ) {
  print("------------------------")
  print((i-1))
  print("------------------------")
  print( cluster.sex.diff[[i]][ cluster.sex.diff[[i]][,5] <= 0.05, ] )
}


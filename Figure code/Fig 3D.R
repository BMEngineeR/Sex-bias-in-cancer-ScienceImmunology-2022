rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)
library(stringr)

# load data

load( file="seurat_data_aggregated_with_annotation.RData" )
all.genes <- rownames(experiment.aggregate)
experiment.aggregate <- ScaleData(experiment.aggregate, features = all.genes)

# effector score

markers <- read.csv( "gene_list_effector.csv" )[,1]

ind <- match( markers, rownames(experiment.aggregate[["RNA"]]@data) )
expr.markers <- t(as.matrix(experiment.aggregate[["RNA"]]@data[ ind[ !is.na(ind) ], ]))
#experiment.aggregate@meta.data$effector <- rowSums(expr.markers)
fit.pca <- princomp(expr.markers, cor=TRUE)
#fit.pca$loadings
experiment.aggregate@meta.data$effector <- fit.pca$scores[,1]

# progenitor exhausted score

markers <- read.csv( "gene_list_prog_exh.csv" )[,1]

ind <- match( markers, rownames(experiment.aggregate[["RNA"]]@data) )
expr.markers <- t(as.matrix(experiment.aggregate[["RNA"]]@data[ ind[ !is.na(ind) ], ]))
#experiment.aggregate@meta.data$prog_exh <- rowSums(expr.markers)
fit.pca <- princomp(expr.markers, cor=TRUE)
#fit.pca$loadings
experiment.aggregate@meta.data$prog_exh <- fit.pca$scores[,1]

# visualization

pdf( "effector.pdf" )
p1 <- FeaturePlot( experiment.aggregate, "effector", col=c("lightyellow","orange","darkred"), min.cutoff="q5", max.cutoff="q90" )
print(p1)
dev.off()

pdf( "prog_exh.pdf" )
p1 <- FeaturePlot( experiment.aggregate, "prog_exh", col=c("lightyellow","orange","darkred"), min.cutoff="q5", max.cutoff="q90" )
print(p1)
dev.off()


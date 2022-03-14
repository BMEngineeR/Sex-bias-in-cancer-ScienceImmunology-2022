rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)
library(RColorBrewer)

# load data

load( file="seurat_data_aggregated_with_annotation.RData" )

# load enrichment analysis results

enrichment.male <- read.csv( file="pvalue_hormone_male.csv",
  header=TRUE, stringsAsFactors=FALSE )
enrichment.female <- read.csv( file="pvalue_hormone_female.csv",
  header=TRUE, stringsAsFactors=FALSE )

# submsaple cells for each of males and females

table(experiment.aggregate@meta.data$sex)
  
experiment.aggregate <- subset( experiment.aggregate, 
  cells=c( 
    base::sample( which( experiment.aggregate@meta.data$sex == "Male" ), 9500, replace=FALSE ),
    base::sample( which( experiment.aggregate@meta.data$sex == "Female" ), 9500, replace=FALSE ) ) )

table(experiment.aggregate@meta.data$sex)

# Visualization of two enrichment analysis p-values on UMAP

loc.male <- which( experiment.aggregate@meta.data$sex == "Male" )
loc.female <- which( experiment.aggregate@meta.data$sex == "Female" )

clust.male <- experiment.aggregate@meta.data$seurat_clusters[ loc.male ]
clust.female <- experiment.aggregate@meta.data$seurat_clusters[ loc.female ]

pos.male <- experiment.aggregate@reductions$umap@cell.embeddings[ loc.male, ]
pos.female <- experiment.aggregate@reductions$umap@cell.embeddings[ loc.female, ]

x.all <- c( pos.male[,1], pos.female[,1] )
y.all <- c( pos.male[,2], pos.female[,2] )

pdf( paste(dir.out,"enrichment_on_UMAP_androgen_interferon.pdf",sep="") )

col.male <- rep( NA, 11 )
for ( i in 1:11 ) {
  col.male[ clust.male == (i-1) ] <- enrichment.male[i,(2+1)]
}
col.female <- rep( NA, 11 )
for ( i in 1:11 ) {
  col.female[ clust.female == (i-1) ] <- enrichment.female[i,(2+1)]
}

col.red.org <- c( col.male, col.female )
col.red <- rep( 1, length(col.red.org) )
for ( gg in 1:4 ) {
  col.red[ col.red.org >= (0.1/4)*(gg-1) & col.red.org <= (0.1/4)*gg ] <- (6-gg)
}

col.male <- rep( NA, 11 )
for ( i in 1:11 ) {
  col.male[ clust.male == (i-1) ] <- enrichment.male[i,(2+4)]
}
col.female <- rep( NA, 11 )
for ( i in 1:11 ) {
  col.female[ clust.female == (i-1) ] <- enrichment.female[i,(2+4)]
}

col.blue.org <- c( col.male, col.female )
col.blue <- rep( 1, length(col.blue.org) )
for ( gg in 1:4 ) {
  col.blue[ col.blue.org >= (0.005/4)*(gg-1) & col.blue.org <= (0.005/4)*gg ] <- (6-gg)
}

col.red[ col.red <= 2 ] <- NA

par( mfrow=c(1,1) )
plot( rbind( pos.male, pos.male, pos.female, pos.female ),
  col=c( brewer.pal(n = 11, name = "RdBu")[ 5 + col.blue[1:9500] ], 
    brewer.pal(n = 11, name = "RdBu")[ 7 - col.red[1:9500] ],
    brewer.pal(n = 11, name = "RdBu")[ 5 + col.blue[9501:19000] ], 
    brewer.pal(n = 11, name = "RdBu")[ 7 - col.red[9501:19000] ] ),
  pch=20, cex=0.1, xlim=c( min(x.all), max(x.all) ), ylim=c( min(y.all), max(y.all) ) )

dev.off()

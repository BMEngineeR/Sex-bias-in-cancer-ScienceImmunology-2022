rm( list=ls() )
#graphics.off()
#library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# options

list.genes <- c( 
  "Xist", "Ddx3y", 
  "Havcr2", "Pdcd1", "Entpd1", "Lag3", "Ctla4", 
  "Cd28", "Cd44", "Cd69", "Tnfrsf4", "Tnfrsf9", "Il2ra", "Icos",
  "Ifng", "Gzmk", "Gzma", "Gzmb", "Prf1",
  "Hif1a", "Id2", "Tcf7", "Slamf6", "Jun", "Fos", "Fosb", "Ikzf2",
  "Sell", "Ccr7", "Csf1", "Ccl3", "Ccl4", "Ccr2", "Ccr5", "Cxcr6", "Bcl2"
  )

# load signature

seurat.marker <- vector( "list", 11 )
for ( i in 1:11 ) {
  marker.i <- read.csv( paste("seurat_list_markers_",(i-1),".csv",sep=""),
    stringsAsFactors=FALSE, header=TRUE )
  seurat.marker[[i]] <- marker.i
}

# load DE genes

DE.clust <- vector( "list", 11 )
for ( i in 1:11 ) {
  print(i)
  DE.i <- read.table( paste("cluster_sex_diff_cluster",(i-1),".txt",sep=""),
    stringsAsFactors=FALSE, header=TRUE )
  DE.clust[[i]] <- DE.i
}

# annotate

col1 <- col2 <- col3 <- col4 <- c()
for ( i in 1:length(DE.clust) ) {
  for ( j in 1:length(list.genes) ) {
    
    loc.ij <- which( seurat.marker[[i]][,1] == list.genes[j] )
    if ( length(loc.ij) >= 1 ) {
      pval <- ( seurat.marker[[i]][loc.ij,3] + seurat.marker[[i]][loc.ij,8] ) / 2
      pval[ pval > 0.5 ] <- 0.5
      pval[ pval < -0.5 ] <- -0.5
      col1 <- c( col1, pval )
    } else {
      col1 <- c( col1, 0 )
    }
    
    loc.ij <- which( rownames(DE.clust[[i]]) == list.genes[j] )
    if ( length(loc.ij) >= 1 ) {
      if ( DE.clust[[i]][loc.ij,2] > 0 ) {
        direction <- 1
      } else {
        direction <- -1
      }
      col2 <- c( col2, direction )
    } else {
      col2 <- c( col2, 0 )
    }
      
    col3 <- c( col3, (i-1) )
    
    col4 <- c( col4, list.genes[j] )
  }
}

col3 <- factor( col3, levels=c(0,5,4,3,8,6,10,1,9,2,7) )
col4 <- factor( col4, levels=list.genes )
pdata <- data.frame( col1, col2, col3, col4 )
colnames(pdata) <- c( "pval", "direction", "cluster", "gene" )

# Fig 3C

dmat <- c()
for ( i in c(7,2,9,1,10,6,8,3,4,5,0) ) {
  dmat <- rbind( dmat, pdata$pval[ pdata$cluster == i ] )
}
rownames(dmat) <- c(7,2,9,1,10,6,8,3,4,5,0) 
colnames(dmat) <- pdata$gene[ pdata$cluster == 0 ]

pheatmap( dmat, color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50), cluster_rows=FALSE, cluster_cols=FALSE, angle_col=45, fontsize_col=20, fontsize_row=20 )

# Fig S2A

dmat <- c()
for ( i in c(7,2,9,1,10,6,8,3,4,5,0) ) {
  dmat <- rbind( dmat, pdata$direction[ pdata$cluster == i ] )
}
rownames(dmat) <- c(7,2,9,1,10,6,8,3,4,5,0) 
colnames(dmat) <- pdata$gene[ pdata$cluster == 0 ]

pheatmap( dmat, color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50), cluster_rows=FALSE, cluster_cols=FALSE, angle_col=45, fontsize_col=20, fontsize_row=20 )

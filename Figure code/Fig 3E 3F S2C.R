rm( list=ls() )
graphics.off()
library(Seurat)
library(cowplot)
library(slingshot)
library(tools)
library(grDevices)
library(RColorBrewer)
library(ggplot2)

# load the Seurat fit

load( file="seurat_data_aggregated_with_annotation.RData" )

# Visualization

p1 <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "sex")
p2 <- DimPlot(experiment.aggregate, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(experiment.aggregate, reduction = "umap", split.by = "sex")

# load into slingshot

umap <- experiment.aggregate@reductions$umap@cell.embeddings
seurat_clusters <- experiment.aggregate@meta.data$seurat_clusters

# select cells

cells.sel <- which( seurat_clusters == 6 |
  seurat_clusters == 10 |
  seurat_clusters == 2 |
  seurat_clusters == 7 |
  seurat_clusters == 1 |
  seurat_clusters == 9 )
umap <- umap[ cells.sel,  ]
seurat_clusters <- seurat_clusters[ cells.sel ]

cells.sel <- which( umap[,1] <= 4 & umap[,2] >=-3 )
umap <- umap[ cells.sel,  ]
seurat_clusters <- seurat_clusters[ cells.sel ]

# run slingshot

fit <- slingshot( data=umap, clusterLabels=seurat_clusters )

pseudotime <- slingPseudotime(fit)

# Fig 3E. clusters

colmap <- rep( "", length(seurat_clusters) )
colmap[ seurat_clusters == 1 ] <- "#db8e00"
colmap[ seurat_clusters == 2 ] <- "#aea200"
colmap[ seurat_clusters == 7 ] <- "#00a6ff"
colmap[ seurat_clusters == 9 ] <- "#ef67eb"
colmap[ seurat_clusters == 6 ] <- "#00bade"
colmap[ seurat_clusters == 10 ] <- "#ff63b6"

plot( umap, col=colmap, pch=18, cex=0.7 )
lines( SlingshotDataSet(fit), lwd=10 )

# Fig 3E. pseudotime

colors <- colorRampPalette(c("lightyellow", "orange", "darkred"))(100)
plotcol <- colors[cut(pseudotime, breaks=100)]

plot( umap, col=plotcol, pch=18, cex=0.7 )
lines( SlingshotDataSet(fit), lwd=10 )

# Fig 3E. marker expression

DefaultAssay(experiment.aggregate) <- "RNA"
scale.data <- as.matrix(GetAssayData(experiment.aggregate, slot = "data"))

list.genes <- c( "Tcf7", "Cd44", "Gzmb", "Prf1", "Havcr2", "Pdcd1" )

pdf( paste(dir.out,"marker_subset.pdf",sep="") )
for ( i in 1:length(list.genes) ) {
  data.i <- scale.data[ rownames(scale.data) == list.genes[i], 
    match( rownames(umap), colnames(scale.data) ) ]
  
  colors <- colorRampPalette(c("lightyellow", "orange", "darkred"))(100)
  plotcol <- colors[cut(data.i, breaks=100)]
  
  plot( umap, col=plotcol, pch=18, cex=0.7, main=list.genes[i] )
  lines( SlingshotDataSet(fit), lwd=10 )
}
dev.off()

# Fig 3F

seurat_clusters <- experiment.aggregate@meta.data$seurat_clusters
sex <- experiment.aggregate@meta.data$sex

sex <- sex[ match( rownames(pseudotime), rownames(experiment.aggregate@meta.data) ) ]
seurat_clusters <- seurat_clusters[ match( rownames(pseudotime), rownames(experiment.aggregate@meta.data) ) ]

pdata <- data.frame( sex, pseudotime )
colnames(pdata) <- c( "sex", "pseudotime" )

p.value <- wilcox.test( pdata$pseudotime ~ pdata$sex )$p.value

p <- ggplot( pdata, aes( x=sex, y=pseudotime, fill=sex ) ) +
  geom_boxplot() + xlab( "Sex" ) + ylab( "Pseudotime" ) +
  ggtitle( paste("Overall (p-value = ",p.value,")",sep="") ) + 
  theme( axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),
         axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
         legend.title=element_text(size=14),legend.text=element_text(size=14),
         panel.background = element_rect(fill = "white", colour="black") )
print(p)

# Fig S2C

data.gzmb <- scale.data[ rownames(scale.data) == "Gzmb", 
  match( rownames(umap), colnames(scale.data) ) ]
data.tcf7 <- scale.data[ rownames(scale.data) == "Tcf7", 
  match( rownames(umap), colnames(scale.data) ) ]

all( names(data.gzmb) == rownames(pseudotime) )

colmap <- rep( "", length(seurat_clusters) )
colmap[ seurat_clusters == 1 ] <- "#db8e00"
colmap[ seurat_clusters == 2 ] <- "#aea200"
colmap[ seurat_clusters == 7 ] <- "#00a6ff"
colmap[ seurat_clusters == 9 ] <- "#ef67eb"
colmap[ seurat_clusters == 6 ] <- "#00bade"
colmap[ seurat_clusters == 10 ] <- "#ff63b6"

plot( pseudotime, data.gzmb, col=colmap, 
  pch=18, cex=0.7, xlab="Pseudotime", ylab="Gzmb" )
lines(loess.smooth( pseudotime[ sex == "Male" ], data.gzmb[ sex == "Male" ], span=4/5 ), lwd=3 )
lines(loess.smooth( pseudotime[ sex == "Female" ], data.gzmb[ sex == "Female" ], span=4/5 ), lwd=3, lty=2 )
legend( "topright", lwd=3, lty=1:2, c( "Male", "Female" ) )

plot( pseudotime, data.tcf7, col=colmap, 
  pch=18, cex=0.7, xlab="Pseudotime", ylab="Tcf7" )
lines(loess.smooth( pseudotime[ sex == "Male" ], data.tcf7[ sex == "Male" ], span=4/5 ), lwd=3 )
lines(loess.smooth( pseudotime[ sex == "Female" ], data.tcf7[ sex == "Female" ], span=4/5 ), lwd=3, lty=2 )
legend( "topleft", lwd=3, lty=1:2, c( "Male", "Female" ) )



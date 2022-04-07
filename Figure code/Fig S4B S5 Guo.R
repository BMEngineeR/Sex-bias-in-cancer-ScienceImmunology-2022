rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)
library(ggplot2)
library(tools)

# load androgen and interferon gene sets

androgen <- read.table( "HALLMARK_ANDROGEN_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
androgen <- toTitleCase(tolower(androgen))
interferon <- read.table( "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
interferon <- toTitleCase(tolower(interferon))

# load data

data.count <- read.table( "GSE99254_NSCLC.TCell.S12346.count.txt", 
  header=TRUE, stringsAsFactors=FALSE )

metadata <- read.table( "cell_annotation.txt", 
  header=TRUE, stringsAsFactors=FALSE, sep="\t" )

# subset for TTC

ind.TTC <- which( metadata$sampleType == "TTC" )
metadata.sub <- metadata[ ind.TTC, ]

# remove filtered cells

ind.filtered <- which( metadata.sub$majorCluster == "filtered" )
metadata.sub <- metadata.sub[ -ind.filtered, ]

# match count data and metadata

ind.sel <- match( gsub( "-", ".", metadata.sub[,1] ), colnames(data.count) )
data.count.sub <- data.count[ , ind.sel ]

data.count.sub <- data.count.sub[ !is.na(data.count$symbol), ]
rownames(data.count.sub) <- data.count$symbol[ !is.na(data.count$symbol) ]

all( gsub( "-", ".", metadata.sub[,1] ) == colnames(data.count.sub) )

# annotate sex

metadata.sub$Sex <- "Female"
metadata.sub$Sex[ metadata.sub$Patient == "P0913" | metadata.sub$Patient == "P0617" | 
    metadata.sub$Patient == "P0616P" | metadata.sub$Patient == "P0619" | 
    metadata.sub$Patient == "P0706" | metadata.sub$Patient == "P1120" ] <- "Male"
rownames(metadata.sub) <- gsub( "-", ".", metadata.sub[,1] )

# construct Seurat dataset

guo <- CreateSeuratObject( counts = data.count.sub, meta.data = metadata.sub, 
  project = "Guo")

# data preprocessing

guo <- NormalizeData( guo, normalization.method = "LogNormalize", scale.factor = 10000 )
guo <- FindVariableFeatures(guo)
guo[["percent.mt"]] <- PercentageFeatureSet( guo, pattern = "^mt-" )
guo <- ScaleData( guo, vars.to.regress = "percent.mt" )

# PCA

guo <- RunPCA( guo, npcs = 100, ndims.print = 1:5, nfeatures.print = 5 )
ElbowPlot( guo, ndims = 100 )
DimHeatmap( guo, dims = c(1:3, 70:75), cells = 500, balanced = TRUE )

# UMAP

guo <- FindNeighbors( guo, reduction = "pca", dims = 1:75, nn.eps = 0.5 )
guo <- FindClusters( guo, resolution = 3, n.start = 10 )
guo <- RunUMAP( guo, dims = 1:75, min.dist = 0.75 )

# select only CD8 cell clusters

guo.CD8 <- guo
guo.CD8 <- subset( guo.CD8, 
  cells=grep( "CD8", guo.CD8@meta.data$majorCluster ) )

# Fig S4B: UMAP

DimPlot(guo.CD8, reduction = "umap", group.by = "majorCluster", split.by = "Sex")
plot_grid(p1)

# Fig S4B: clustering results

table( guo@meta.data$sex )
table( guo@meta.data$majorCluster )
clust.M <- guo@meta.data$majorCluster[ guo@meta.data$Sex == "Male" ]
clust.F <- guo@meta.data$majorCluster[ guo@meta.data$Sex == "Female" ]
cluster.uniq <- unique(guo@meta.data$majorCluster)
freq.table <- matrix( NA, 2, length(cluster.uniq) )
for ( i in 1:length(cluster.uniq) ) {
  freq.table[1,i] <- round( 100 * length(which( clust.M == cluster.uniq[i] )) / length(clust.M) ) / 100
  freq.table[2,i] <- round( 100 * length(which( clust.F == cluster.uniq[i] )) / length(clust.F) ) / 100
}
rownames(freq.table) <- c( "Male", "Female" )
colnames(freq.table) <- cluster.uniq
freq.table

freq.table2 <- freq.table[,1:7]
freq.table2[1,] <- freq.table2[1,] / sum(freq.table2[1,])
freq.table2[2,] <- freq.table2[2,] / sum(freq.table2[2,])
round( 100 * freq.table2 ) / 100

# extract only the cells corresponding to male

guo.male <- guo
guo.male <- subset( guo.male, 
  cells=which( guo.male@meta.data$Sex == "Male" ) )

# correlation analysis

raw.data <- as.matrix(GetAssayData(guo.male, slot = "data"))
loc.androgen <- match( toupper(androgen), rownames(raw.data) )
data.androgen <- unlist(raw.data[ loc.androgen[ !is.na(loc.androgen) ], ])
loc.interferon <- match( toupper(interferon), rownames(raw.data) )
data.interferon <- unlist(raw.data[ loc.interferon[ !is.na(loc.interferon) ], ])

cor.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
gene.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
for ( i in 1:nrow(data.androgen) ) {
  for ( j in 1:nrow(data.interferon) ) {
    x <- data.androgen[i,]
    y <- data.interferon[j,]
    x1 <- x[ -which( x == 0 & y == 0 ) ]
    y1 <- y[ -which( x == 0 & y == 0 ) ]
    cor.mat[i,j] <- cor( x1, y1 )
    gene.mat[i,j] <- paste( rownames(data.androgen)[i], rownames(data.interferon)[j] )
  }
}
boxplot(as.vector(cor.mat))
summary(as.vector(cor.mat))

data.male <- data.frame( as.vector(gene.mat), as.vector(cor.mat) )
data.male <- data.male[ order(data.male[,2]), ]

cors <- list()
cors[["Male"]] <- as.vector(cor.mat)

# extract only the cells corresponding to female

guo.female <- guo
guo.female <- subset( guo.female, 
  cells=which( guo.female@meta.data$Sex == "Female" ) )

# correlation analysis

raw.data <- as.matrix(GetAssayData(guo.female, slot = "data"))
loc.androgen <- match( toupper(androgen), rownames(raw.data) )
data.androgen <- unlist(raw.data[ loc.androgen[ !is.na(loc.androgen) ], ])
loc.interferon <- match( toupper(interferon), rownames(raw.data) )
data.interferon <- unlist(raw.data[ loc.interferon[ !is.na(loc.interferon) ], ])

cor.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
gene.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
for ( i in 1:nrow(data.androgen) ) {
  for ( j in 1:nrow(data.interferon) ) {
    x <- data.androgen[i,]
    y <- data.interferon[j,]
    x1 <- x[ -which( x == 0 & y == 0 ) ]
    y1 <- y[ -which( x == 0 & y == 0 ) ]
    cor.mat[i,j] <- cor( x1, y1 )
    gene.mat[i,j] <- paste( rownames(data.androgen)[i], rownames(data.interferon)[j] )
  }
}
boxplot(as.vector(cor.mat))
summary(as.vector(cor.mat))

data.female <- data.frame( as.vector(gene.mat), as.vector(cor.mat) )
data.female <- data.female[ order(data.female[,2]), ]

cors[["Female"]] <- as.vector(cor.mat)

# Fig S5B

cors[["Male"]][ cors[["Male"]] > 0.5 ] <- NA
cors[["Female"]][ cors[["Female"]] > 0.5 ] <- NA

pdata <- data.frame(
  c( rep( "Male", length(cors[["Male"]]) ), rep( "Female", length(cors[["Female"]]) ) ),
  c( cors[["Male"]], cors[["Female"]] ) 
)
colnames(pdata) <- c( "Sex", "Corr" )

p <- ggplot( pdata, aes( x=Sex, y=Corr, fill=Sex ) ) + 
  geom_boxplot() + ylab( "Correlation coefficient" ) + ylim( -1, 1 ) +
  geom_hline( yintercept=0, linetype="dashed", color="red", size=2 ) +
  theme(axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
        legend.title=element_text(size=14),legend.text=element_text(size=14))
print(p)


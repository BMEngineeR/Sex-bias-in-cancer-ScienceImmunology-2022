rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)
library(ggplot2)

# load androgen and interferon gene sets

androgen <- read.table( "HALLMARK_ANDROGEN_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
androgen <- toupper(androgen)
interferon <- read.table( "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
interferon <- toupper(interferon)

# load data

data.count <- read.table( "GSE123813_bcc_scRNA_counts.txt", 
  header=TRUE, stringsAsFactors=FALSE )

metadata.all <- read.table( "bcc_all_metadata.txt", 
  stringsAsFactors=FALSE, sep="\t" )
metadata.tcell <- read.table( "bcc_tcell_metadata.txt", 
  stringsAsFactors=FALSE, sep="\t" )

# subset for t-cell, pre-treatment, & CD8

ind.pre <- which( metadata.tcell[,2] == "pre")
ind.CD8 <- grep( "CD8", metadata.tcell[,3] )

metadata.tcell.sub <- metadata.tcell[ intersect( ind.pre, ind.CD8 ), ]

ind.sel <- match( rownames(metadata.tcell.sub), colnames(data.count) )
data.count.sub <- data.count[ , ind.sel ]

# construct Seurat dataset

yost <- CreateSeuratObject( counts = data.count.sub, meta.data = metadata.tcell.sub, 
  project = "Yost")

# data preprocessing

yost <- NormalizeData( yost, normalization.method = "LogNormalize", scale.factor = 10000 )
yost <- FindVariableFeatures(yost)
yost[["percent.mt"]] <- PercentageFeatureSet( yost, pattern = "^mt-" )
yost <- ScaleData( yost, vars.to.regress = "percent.mt" )

# PCA

yost <- RunPCA( yost, npcs = 100, ndims.print = 1:5, nfeatures.print = 5 )
ElbowPlot( yost, ndims = 100 )
DimHeatmap( yost, dims = c(1:3, 70:75), cells = 500, balanced = TRUE )

# UMAP

yost <- FindNeighbors( yost, reduction = "pca", dims = 1:75, nn.eps = 0.5 )
yost <- FindClusters( yost, resolution = 3, n.start = 10 )
yost <- RunUMAP( yost, dims = 1:75, min.dist = 0.75 )

# heatmap for sex genes

DefaultAssay(yost) <- "RNA"
Idents(yost) <- "patient"
DoHeatmap(yost, features = c("XIST", "DDX3Y"), slot="scale.data" )

# annotate for sex

patient.sex <- rep( "Male", length(yost@meta.data$patient) )
patient.sex[ yost@meta.data$patient == "su006" | yost@meta.data$patient == "su007" |
    yost@meta.data$patient == "su012" ] <- "Female"
yost$sex <- patient.sex
Idents(yost) <- "sex"

# Fig S4A: UMAP

DimPlot(yost, reduction = "umap", group.by = "cluster", split.by = "sex")
plot_grid(p1)

# Fig S4A: clustering results

table( yost@meta.data$sex )
table( yost@meta.data$cluster )
clust.M <- yost@meta.data$cluster[ yost@meta.data$sex == "Male" ]
clust.F <- yost@meta.data$cluster[ yost@meta.data$sex == "Female" ]
cluster.uniq <- unique(yost@meta.data$cluster)
freq.table <- matrix( NA, 2, length(cluster.uniq) )
for ( i in 1:length(cluster.uniq) ) {
  freq.table[1,i] <- round( 100 * length(which( clust.M == cluster.uniq[i] )) / length(clust.M) ) / 100
  freq.table[2,i] <- round( 100 * length(which( clust.F == cluster.uniq[i] )) / length(clust.F) ) / 100
}
rownames(freq.table) <- c( "Male", "Female" )
colnames(freq.table) <- cluster.uniq
freq.table

# extract only the cells corresponding to male

yost.male <- yost
yost.male <- subset( yost.male, 
  cells=which( yost.male@meta.data$sex == "Male" ) )

DefaultAssay(yost.male) <- "RNA"
Idents(yost.male) <- "cluster"
all.genes <- rownames(yost.male)
yost.male <- ScaleData(yost.male, features = all.genes)

# correlation analysis

raw.data <- as.matrix(GetAssayData(yost.male, slot = "data"))
loc.androgen <- match( androgen, rownames(raw.data) )
data.androgen <- unlist(raw.data[ loc.androgen[ !is.na(loc.androgen) ], ])
loc.interferon <- match( interferon, rownames(raw.data) )
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

yost.female <- yost
yost.female <- subset( yost.female, 
  cells=which( yost.female@meta.data$sex == "Female" ) )

DefaultAssay(yost.female) <- "RNA"
Idents(yost.female) <- "cluster"
all.genes <- rownames(yost.female)
yost.female <- ScaleData(yost.female, features = all.genes)
dev.off()

# correlation analysis

raw.data <- as.matrix(GetAssayData(yost.female, slot = "data"))
loc.androgen <- match( androgen, rownames(raw.data) )
data.androgen <- unlist(raw.data[ loc.androgen[ !is.na(loc.androgen) ], ])
loc.interferon <- match( interferon, rownames(raw.data) )
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

pdata <- data.frame(
  c( rep( "Male", length(cors[["Male"]]) ), rep( "Female", length(cors[["Female"]]) ) ),
  c( cors[["Male"]], cors[["Female"]] )
)
colnames(pdata) <- c( "Sex", "Correlation" )

pdata <- pdata[ !is.na(pdata[,2]), ]
pdata <- pdata[ pdata[,2] < 0.5, ]

pdf( paste(dir.out,"corr_androgen_interferon_ggplot2.pdf",sep="") )
p <- ggplot( pdata, aes( x=Sex, y=Correlation, fill=Sex ) ) +
  geom_boxplot() + xlab( "Sex" ) + ylab( "Correlation coefficients" ) + ylim( -1, 1 ) +
  geom_hline( yintercept=0, linetype="dashed", color="red", size=2 ) +
  theme( axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),
         axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
         legend.title=element_text(size=14),legend.text=element_text(size=14),
         panel.background = element_rect(fill = "white", colour="black") )
print(p)
dev.off()


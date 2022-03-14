rm( list=ls() )
graphics.off()
library(monocle)
library(tools)
library(ggplot2)

# load androgen and interferon gene sets

androgen <- read.table( "HALLMARK_ANDROGEN_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
androgen <- toTitleCase(tolower(androgen))

interferon <- read.table( "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
interferon <- toTitleCase(tolower(interferon))

# load data

count1 <- read.table( "GSE131535_matrix.mtx", 
  header=FALSE, stringsAsFactors=FALSE, skip=3 )

gene1 <- read.table( "GSE131535_genes.tsv", 
  header=FALSE, stringsAsFactors=FALSE )

cell1 <- read.table( "GSE131535_barcodes.tsv", 
  header=FALSE, stringsAsFactors=FALSE, sep="\t" )

# processing

data.count1 <- matrix( 0, nrow(gene1), nrow(cell1) )
rownames(data.count1) <- gene1[,1]
colnames(data.count1) <- cell1[,1]
for ( i in 1:nrow(count1) ) {
  data.count1[ count1[i,1], count1[i,2] ] <- count1[i,3]
}

# load into monocle

expr.matrix <- data.count1
cell.metadata <- cell1

gene.annotation <- data.frame(gene1)
rownames(gene.annotation) <- gene1[,1]
colnames(gene.annotation) <- c( "gene_short_name", "symbol" )

rownames(cell.metadata) <- cell.metadata[,1]

# select genes with large variation

sd.log.count <- apply( log10(expr.matrix+1), 1, sd )
loc.sel <- c( 1:nrow(expr.matrix) )[ order( sd.log.count, decreasing=TRUE ) ][1:2000]

expr.matrix.org <- expr.matrix

expr.matrix <- expr.matrix[ loc.sel, ]
gene.annotation <- gene.annotation[ loc.sel, ]

# Choosing a distribution and construct an object

pd <- new("AnnotatedDataFrame", data = cell.metadata)
fd <- new("AnnotatedDataFrame", data = gene.annotation)

cd8t <- newCellDataSet(expr.matrix,
    phenoData = pd,
    featureData = fd,
    expressionFamily=negbinomial.size())

# Estimate size factors and dispersions

cd8t <- estimateSizeFactors(cd8t)
cd8t <- estimateDispersions(cd8t)

# cell classification

disp_table <- dispersionTable(cd8t)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cd8t <- setOrderingFilter(cd8t, unsup_clustering_genes$gene_id)
plot_ordering_genes(cd8t)

plot_pc_variance_explained(cd8t, return_all = F)

cd8t <- reduceDimension(cd8t, max_components = 2, num_dim = 6,
  reduction_method = 'tSNE', verbose = T)
cd8t <- clusterCells(cd8t, num_clusters = 6)

plot_cell_clusters(cd8t)

# Constructing Single Cell Trajectories: (1) order genes

clustering_DEG_genes <-
    differentialGeneTest(cd8t,
      fullModelFormulaStr = '~Cluster',
      cores = 10)

cd8t_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

cd8t <-setOrderingFilter(cd8t,
    ordering_genes = cd8t_ordering_genes)

# Constructing Single Cell Trajectories: (2) reduce data dimensionality 

cd8t <- reduceDimension(cd8t, max_components = 2,
    method = 'DDRTree')

# Constructing Single Cell Trajectories: (3) order cells along the trajectory 

cd8t <- orderCells(cd8t)

# Fig 6C

androgen.keys <- c( "Dbi", "Myl12a", "Sat1", "Tmem50a" )

loc.sel <- match( androgen.keys, gene1[,2] )
androgen.sel <- androgen[ !is.na(loc.sel) ]
loc.sel <- loc.sel[ !is.na(loc.sel) ]
ensembl.id <- gene1[ loc.sel, 1 ]
expr.sel <- expr.matrix.org[ ensembl.id, ]

data.state.effector <- expr.sel[ , pData(cd8t)$State == 3 ] + 1
data.state.exausted <- expr.sel[ , pData(cd8t)$State == 1 ] + 1

pdata <- c()
for ( i in 1:nrow(data.state.effector) ) {
  pdata.i <- data.frame( rep( androgen.keys[i], length(data.state.effector[i,]) + length(data.state.exausted[i,]) ), 
    c( rep( "Eff", length(data.state.effector[i,]) ), rep( "Ex", length(data.state.exausted[i,]) ) ),
    c( data.state.effector[i,], data.state.exausted[i,] ) )
  names(pdata.i) <- c( "Gene", "Cell", "Expr" )
  pdata <- rbind( pdata, pdata.i )
}

p <- ggplot( pdata, aes( x=Cell, y=Expr, fill=Cell ) ) + facet_grid( . ~ Gene ) +
  scale_y_continuous( trans='log10' ) +
  geom_boxplot() + xlab( "Cell branch" ) + ylab( "Gene expression" ) +
  #geom_hline( yintercept=0, linetype="dashed", color="red", size=2 ) +
  theme( axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),
		axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
		legend.title=element_text(size=14),legend.text=element_text(size=14),
		strip.background=element_rect( colour="black", fill="white", size=1 ),
		strip.text.x=element_text( size=20 ),
		panel.background = element_rect(fill = "white", colour="black") )
print(p)

# Fig 6B

loc.sel <- match( androgen, gene1[,2] )
androgen.sel <- androgen[ !is.na(loc.sel) ]
loc.sel <- loc.sel[ !is.na(loc.sel) ]
ensembl.id <- gene1[ loc.sel, 1 ]
expr.sel <- expr.matrix.org[ ensembl.id, ]

data.state.branch1 <- log10( expr.sel[ , pData(cd8t)$State == 1 ] + 1 )
data.state.branch2 <- log10( expr.sel[ , pData(cd8t)$State == 2 ] + 1 )
data.state.branch3 <- log10( expr.sel[ , pData(cd8t)$State == 3 ] + 1 )

mat.ratio <- matrix( NA, nrow(expr.sel), length(pData(cd8t)$State) )
for ( i in 1:nrow(expr.sel) ) {
  mat.ratio[ i, pData(cd8t)$State == 2 ] <- 0
  mat.ratio[ i, pData(cd8t)$State == 1 ] <- mean(data.state.branch1[i,]) - mean(data.state.branch2[i,])
  mat.ratio[ i, pData(cd8t)$State == 3 ] <- mean(data.state.branch3[i,]) - mean(data.state.branch2[i,])
}
mat.ratio[ is.nan(mat.ratio) ] <- 0

pData(cd8t)$Ratio <- apply( mat.ratio, 2, mean )

plot_cell_trajectory(cd8t, color_by = "Ratio")

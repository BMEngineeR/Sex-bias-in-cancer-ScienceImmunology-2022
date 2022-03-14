rm( list=ls() )
graphics.off()
library(reshape2)
library(tools)

# options

list.geneset <- c( "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE" )

# load signatures

gene.hormone <- vector( "list", length(list.geneset) )
for ( i in 1:length(list.geneset) ) {
  list.i <- read.table( paste(list.geneset[i],".txt",sep=""),
      stringsAsFactors=FALSE, header=FALSE, skip=2 )
  gene.hormone[[i]] <- list.i[,1]
}

# load gene cluster markers

seurat.marker.male <- vector( "list", 11 )
for ( i in 1:11 ) {
  print(i)
  marker.i <- read.csv( paste("seurat_list_markers_",(i-1),"_male.csv",sep=""),
    stringsAsFactors=FALSE, header=TRUE )
  seurat.marker.male[[i]] <- marker.i
}

seurat.marker.female <- vector( "list", 11 )
for ( i in 1:11 ) {
  print(i)
  marker.i <- read.csv( paste("seurat_list_markers_",(i-1),"_female.csv",sep=""),
    stringsAsFactors=FALSE, header=TRUE )
  seurat.marker.female[[i]] <- marker.i
}

# enrichment analysis

p.hormone.male <- matrix( NA, 11, length(list.geneset) )
for ( j in 1:length(list.geneset) ) {
  for ( i in 1:11 ) {
    gene1 <- gene.hormone[[j]]
    gene2 <- toupper(seurat.marker.male[[i]][,1])
    
    n.ti <- length(intersect( gene1, gene2 ))
    print(i)
    print( seurat.marker.male[[i]][ match( intersect( gene1, gene2 ), toupper(seurat.marker.male[[i]][,1]) ), ] )
    n.t <- length(gene1)
    n.i <- length(gene2)
    n.all <- 10000
    
    p.hormone.male[i,j] <- 
      0.5 * dhyper( n.ti, n.t, n.all - n.t, n.i ) +
      phyper( n.ti, n.t, n.all - n.t, n.i, lower.tail=FALSE )
  }
}

p.hormone.female <- matrix( NA, 11, length(list.geneset) )
for ( j in 1:length(list.geneset) ) {
  for ( i in 1:11 ) {
    gene1 <- gene.hormone[[j]]
    gene2 <- toupper(seurat.marker.female[[i]][,1])
    
    n.ti <- length(intersect( gene1, gene2 ))
    print(i)
    print( seurat.marker.female[[i]][ match( intersect( gene1, gene2 ), toupper(seurat.marker.female[[i]][,1]) ), ] )
    n.t <- length(gene1)
    n.i <- length(gene2)
    n.all <- 10000
    
    p.hormone.female[i,j] <- 
      0.5 * dhyper( n.ti, n.t, n.all - n.t, n.i ) +
      phyper( n.ti, n.t, n.all - n.t, n.i, lower.tail=FALSE )
  }
}

# summary

seurat.ptable.male <- data.frame( paste("cluster",0:10,sep=""), p.hormone.male )
colnames(seurat.ptable.male) <- c( "Seurat cluster", paste("p-value for ",list.geneset,sep="") )
seurat.ptable.male

seurat.ptable.female <- data.frame( paste("cluster",0:10,sep=""), p.hormone.female )
colnames(seurat.ptable.female) <- c( "Seurat cluster", paste("p-value for ",list.geneset,sep="") )
seurat.ptable.female

write.csv( seurat.ptable.male, file=paste(dir.out,"pvalue_hormone_male.csv",sep="") )

write.csv( seurat.ptable.female, file=paste(dir.out,"pvalue_hormone_female.csv",sep="") )


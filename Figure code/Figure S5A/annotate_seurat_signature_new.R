rm( list=ls() )
graphics.off()
library(reshape2)
library(tools)

# options

list.geneset <- c( "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_INTERFERON_ALPHA_RESPONSE" )

#dir.sig <- "/Users/dongjunchung/Dropbox/USB/research - MC23. Zihai Li/6. single cell RNA-seq/4. downstream analysis DE/KEGG ribosome/"
#dir.seurat <- "/Users/dongjunchung/Dropbox/USB/research - MC23. Zihai Li/6. single cell RNA-seq/3. Seurat_subset/results_new_DE/"
#dir.out <- "/Users/dongjunchung/Dropbox/USB/research - MC23. Zihai Li/6. single cell RNA-seq/4. downstream analysis DE/"

dir.sig <- "C:/Users/dchung/Dropbox/USB/research - MC23. Zihai Li/6. single cell RNA-seq/4. downstream analysis DE/hormone/"
dir.seurat <- "C:/Users/dchung/Dropbox/USB/research - MC23. Zihai Li/6. single cell RNA-seq/3. Seurat_subset/results_new/"
dir.out <- "C:/Users/dchung/Dropbox/USB/research - MC23. Zihai Li/6. single cell RNA-seq/4. downstream analysis DE/"

# load signatures

gene.hormone <- vector( "list", length(list.geneset) )
for ( i in 1:length(list.geneset) ) {
  list.i <- read.table( paste(dir.sig,list.geneset[i],".txt",sep=""),
      stringsAsFactors=FALSE, header=FALSE, skip=2 )
  gene.hormone[[i]] <- list.i[,1]
}

# load gene cluster markers

seurat.marker <- vector( "list", 11 )
for ( i in 1:11 ) {
  print(i)
  marker.i <- read.csv( paste(dir.seurat,"seurat_list_markers_",(i-1),".csv",sep=""),
    stringsAsFactors=FALSE, header=TRUE )
  seurat.marker[[i]] <- marker.i
}

# enrichment analysis

p.hormone <- matrix( NA, 11, length(list.geneset) )
for ( j in 1:length(list.geneset) ) {
  for ( i in 1:11 ) {
    gene1 <- gene.hormone[[j]]
    gene2 <- toupper(seurat.marker[[i]][,1])
    
    n.ti <- length(intersect( gene1, gene2 ))
    print(n.ti)
    n.t <- length(gene1)
    n.i <- length(gene2)
    n.all <- 10000
    
    p.hormone[i,j] <- 
      0.5 * dhyper( n.ti, n.t, n.all - n.t, n.i ) +
      phyper( n.ti, n.t, n.all - n.t, n.i, lower.tail=FALSE )
  }
}

# summary

seurat.ptable <- data.frame( paste("cluster",0:10,sep=""), p.hormone )
colnames(seurat.ptable) <- c( "Seurat cluster", paste("p-value for ",list.geneset,sep="") )
seurat.ptable

#write.csv( seurat.ptable, file=paste(dir.out,"pvalue_hormone_all.csv",sep="") )

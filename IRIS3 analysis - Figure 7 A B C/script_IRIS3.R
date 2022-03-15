library(Seurat)
# load Seurat object from previous script named experiment.aggregate
# alternatively, the file can be loaded from here https://cloud.osubmi.com/downloadFiles/seurat_data_aggregated_with_annotation.RData
# and then load the Rdata.
# load("seurat_data_aggregated_with_annotation.RData")
male.object <- subset(experiment.aggregate,subset = sex =="Male")
female.object <- subset(experiment.aggregate,subset = sex =="Female")
# shared  SP3 SP5 SP2 SP1 WT1 ELF5 KLF15 
# male ZNF128 SMCA5 ZNF148 MAZ SALL1 ZBT17 
# female AP2B E2F3 E2F6 SALL4
DefaultAssay(experiment.aggregate) <- "RNA"


# generate 
female.matrix <- as.data.frame(cbind(ID= rownames(female.object),
                                     as.matrix(female.object@assays$RNA@data))
)
female.label <- as.data.frame(cbind(ID =names(female.object$seurat_clusters),
                                    cluster = as.character(female.object$seurat_clusters)
))
identical(colnames(female.matrix)[-1],as.character(female.label$ID))
# 
male.matrix <- as.data.frame(cbind(ID= rownames(male.object),
                                   as.matrix(male.object@assays$RNA@data))
)
male.label <- as.data.frame(cbind(ID =names(male.object$seurat_clusters),
                                  cluster = as.character(male.object$seurat_clusters)
))

identical(colnames(male.matrix)[-1],as.character(male.label$ID))

write.table(female.matrix,file = "female.matrix.txt",quote = F,row.names = F,sep = "\t")
write.table(female.label,file= "female.label.txt",quote = F,row.names = F,sep = "\t")
write.table(male.label,file= "male.label.txt",quote = F,row.names = F,sep = "\t")
write.table(male.matrix,file = "male.matrix.txt",quote = F,row.names = F,sep = "\t")
# after creating these files, upload to IRIS3.


# shared TF by male and female
upstream_10.TFcodinggene <- c("Sp3","Sp2","Sp1","Zfp281","Smarca5","Zfp148","Maz","Zbtb17","E2f3","E2f6","Tcf7")
# male correlation
my.object <- male.object
my.object
my.table <- c()
for (i in 1:length(upstream_10.TFcodinggene)){
  tmp.name <- upstream_10.TFcodinggene[i]
  tcf7._expression.cell <- colnames(my.object)[my.object@assays$RNA@data["Tcf7",]>0]
  TFs._expression.cell <- colnames(my.object)[my.object@assays$RNA@data[tmp.name,]>0]
  overlap.cells <- intersect(tcf7._expression.cell,TFs._expression.cell)
  mean.expression.total <- mean(my.object@assays$RNA@counts[tmp.name,tcf7._expression.cell ])
  tcf7.name <- "Tcf7"
  my.result.all <- cor.test(my.object@assays$RNA@data[tmp.name,overlap.cells ],my.object@assays$RNA@data[tcf7.name,overlap.cells],method = "pearson")
  all.cor <- my.result.all$estimate
  all.p <- my.result.all$p.value
  tmp.table <- c(tmp.name,
                 mean.expression.total,
                 all.cor,
                 all.p
  )
  my.table <- rbind(my.table,tmp.table)
}
my.table <- cbind(my.table,c(as.numeric(as.character(my.table[,4])) * nrow(my.table)))
my.table[,5] <- ifelse(as.numeric(as.character(my.table[,5]))>1, 1,my.table[,5])
colnames(my.table) <- c("TF","Avg.exp",
                        "correlation","p-value","p-adj")

# male and female by different cluster
# use "Sp3","Sp2","Sp1","Zfp281","Smarca5","Zfp148","Maz","Zbtb17","E2f3","E2f6","Tcf7"
# use male.object and female.object
my.object <- male.object
tmp.name <- "Sp3"
my.new.table <- c()
my.cluster <- c(1,2,7,9)
for (i in 1:4){
  cell.name <- names(my.object$seurat_clusters)[my.object$seurat_clusters==my.cluster[i]]
  Idents(my.object) <- my.object$seurat_clusters
  my.tmp.object <- subset(my.object,idents = my.cluster[i])
  mean.expression.total <- mean(my.tmp.object@assays$RNA@data[tmp.name,])
  tcf7.name <- "Tcf7"
  my.result.all <- cor.test(my.tmp.object@assays$RNA@data[tmp.name,],my.tmp.object@assays$RNA@data[tcf7.name,],method = "spearman")
  all.cor <- my.result.all$estimate
  all.p <- my.result.all$p.value
  tmp.table <- c(paste0(c("cluster",my.cluster[i]),collapse = "_"),
                 mean.expression.total,
                 all.cor,
                 all.p
  )
  my.new.table<- rbind(my.new.table,tmp.table)
}

colnames(my.new.table) <- c("cluster","Avg.exp",
                            "correlation","p-value")


# scatter plot 
# change object name by male.object female.object
my.object <- female.object
my.cluster <- c(1,2,7,9)
for (i in 1:4){
  for (j in 1:10){
    cell.name <- names(my.object$seurat_clusters)[my.object$seurat_clusters==my.cluster[i]]
    Idents(my.object) <- my.object$seurat_clusters
    my.tmp.object <- subset(my.object,idents = my.cluster[i])
    tmp.name <- upstream_10.TFcodinggene[j]
    tcf7.name <- "Tcf7"
    
    tcf7._expression.cell <- colnames(my.tmp.object)[my.tmp.object@assays$RNA@data["Tcf7",]>0]
    TFs._expression.cell <- colnames(my.tmp.object)[my.tmp.object@assays$RNA@data[tmp.name,]>0]
    overlap.cells <- intersect(tcf7._expression.cell,TFs._expression.cell)
    
    y.cor<- my.tmp.object@assays$RNA@data[tmp.name,overlap.cells]
    x.cor <- my.tmp.object@assays$RNA@data[tcf7.name,overlap.cells]
    jpeg(paste0("cluster_",i,"_scatter_plot_female_",tmp.name,"_",tcf7.name,".jpeg"))
    plot(x.cor,y.cor,pch = 20,main = paste0("cluster: ",i," scatter plot (female)"),ylab = tmp.name,xlab =tcf7.name )
    abline(lm(x.cor~y.cor), col="red")
    dev.off()
  }
}



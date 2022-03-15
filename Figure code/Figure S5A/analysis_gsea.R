# Implement gene set enrichment analysis
# Dr. Sundi's scRNA-seq data
library(tidyverse)
library(Seurat)

# Interested in the androgen response gene set sent by Tong
ar_set <- read.table("data/HALLMARK_ANDROGEN_RESPONSE.txt") %>%
  pull(V1)

# Load CD8T Seurat object 
load("data/integrated_processed_CD8T.RData")

umap_CD8 <- data.frame(UM1 = exp_agg_CD8@reductions$umap@cell.embeddings[,1],
                       UM2 = exp_agg_CD8@reductions$umap@cell.embeddings[,2],
                       Cluster = as.numeric(exp_agg_CD8$seurat_clusters) - 1)

# Find markers for each cluster
Idents(exp_agg_CD8) <- "seurat_clusters"
cd8t_markers <- FindAllMarkers(exp_agg_CD8)
write.csv(cd8t_markers,
          file = "data/CD8T_cluster_markers.csv",
          quote = FALSE,
          row.names = FALSE)

# Conduct enrichment analysis
## Cluster p-value data frame
n_clust <- length(unique(exp_agg_CD8$seurat_clusters)) # number of clusters
gsea_df <- data.frame(
  Cluster = 0:(n_clust-1),
  Pvalue = rep(0,n_clust)
)

# loop through clusters
for(i in 1:n_clust)
{
  # set of genes for cluster i-1
  genes_i <- cd8t_markers %>%
    filter(cluster == i-1) %>%
    pull(gene) %>%
    toupper()
  
  # number of intersections with AR set
  n.ti <- length(intersect(ar_set, genes_i))
  print(n.ti)
  n.t <- length(ar_set)
  n.i <- length(genes_i)
  n.all <- 10000
  
  # hyper geometric p-value
  gsea_df$Pvalue[i] <- 
    0.5 * dhyper( n.ti, n.t, n.all - n.t, n.i ) +
    phyper( n.ti, n.t, n.all - n.t, n.i, lower.tail=FALSE )
}

# compute -log10 pvalues
gsea_df <- gsea_df %>%
  mutate(ml10_Pvalue = -log10(Pvalue))

write.csv(gsea_df,
          file = "data/AR_CD8T_GSEA.csv",
          quote = FALSE,
          row.names = FALSE)

umap_CD8 <- inner_join(umap_CD8,gsea_df,by = "Cluster")

enrich_plot_bar <- ggplot(data = gsea_df, aes(x = Cluster, y = ml10_Pvalue, fill = ml10_Pvalue)) + 
  geom_col() + 
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_y_continuous(expand = c(0,0)) + 
  ggtitle("AR Enrichment in CD8T")
ggsave(enrich_plot_bar,
       filename = "figures/CD8T_enrichment_bar_plot.pdf",
       height = 6,
       width = 8)

enrich_plot <- ggplot(data = umap_CD8, aes(x = UM1, y = UM2, color = ml10_Pvalue)) + 
  geom_point(size = 2.5) + 
  theme_classic()+ 
  scale_color_continuous(name = "-log10(Pvalue)")+ 
  ggtitle("AR Enrichment in CD8T")
ggsave(enrich_plot,
       filename = "figures/CD8T_enrichment_plot_dot.pdf",
       height = 6,
       width = 8)

enrich_plot_hex <- ggplot(data = umap_CD8, aes(x = UM1, y = UM2, z = ml10_Pvalue)) + 
  stat_summary_hex() + 
  theme_classic() + 
  scale_fill_continuous(name = "-log10(Pvalue)")+ 
  ggtitle("AR Enrichment in CD8T")
ggsave(enrich_plot_hex,
       filename = "figures/CD8T_enrichment_bar_hex.pdf",
       height = 6,
       width = 8)






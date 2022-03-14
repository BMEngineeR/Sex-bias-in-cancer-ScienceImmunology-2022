# load packages
library(ggplot2)
library(ggrepel)
library(dplyr)
library(Seurat)

# load pre-integrated Seurat object
# load("data/Murine_BBN/integrated_seurat.RData")
# DefaultAssay(exp_agg) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
# exp_agg <- ScaleData(exp_agg, verbose = FALSE)
# exp_agg <- RunPCA(exp_agg, npcs = 30, verbose = FALSE)
# 
# # umap and Clustering
# exp_agg <- RunUMAP(exp_agg, reduction = "pca", dims = 1:20)
# exp_agg <- FindNeighbors(exp_agg, reduction = "pca", dims = 1:20)
# exp_agg <- FindClusters(exp_agg, resolution = 0.5)
# 
# # save processed seurat object
# save(exp_agg,
#      file = "data/Murine_BBN/integrated_processed_seurat.RData")
# 
# # save counts for singleR
# counts <- exp_agg[["RNA"]]@counts
# metadata <- exp_agg@meta.data
# save(counts,metadata,file = "data/Murine_BBN/counts.RData")

load("data/Murine_BBN/integrated_processed_seurat.RData")
exp_agg[["percent.mt"]] <- PercentageFeatureSet(exp_agg, pattern = "^Mt-", assay = "RNA")
exp_agg$samp <- "integrated"
Idents(exp_agg) <- "samp"
qc_plot <- VlnPlot(exp_agg, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0) 
ggsave(qc_plot,
       filename = "figures/qc_plot_all.pdf",
       height = 6,
       width = 8)
qc_plot_pt <- VlnPlot(exp_agg, 
                      features = c("nFeature_RNA", "nCount_RNA"), 
                      ncol = 3, pt.size = 0.05) 
ggsave(qc_plot_pt,
       filename = "figures/qc_plot_pt_all.pdf",
       height = 6,
       width = 8)

exp_agg$group_week <- paste(exp_agg$group,exp_agg$week,sep = "_")

# generate umap plot with clusters
p1 <- Seurat::DimPlot(exp_agg, reduction = "umap", label = TRUE) + 
  theme(legend.position = "none")
ggsave(filename = "figures/umap_clusters.pdf",
       plot = p1,
       device = "pdf",
       width = 8,
       height = 7)

p2 <- DimPlot(exp_agg, reduction = "umap", label = TRUE, split.by = "group") + 
  theme(legend.position = "none")
ggsave(filename = "figures/umap_BBNvsControl.pdf", 
       plot = p2,
       device = "pdf",
       width = 8,
       height = 7)

p3 <- DimPlot(exp_agg, reduction = "umap", label = TRUE, split.by = "group_week") + 
  theme(legend.position = "none")
ggsave(filename = "figures/umap_BBNvsControl_week.pdf", 
       plot = p3,
       device = "pdf",
       width = 8,
       height = 7)

# annotations <- read.csv("figures/cell_type_cluster_annotations.csv") %>%
#   dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))
# exp_agg@meta.data <- dplyr::inner_join(exp_agg@meta.data,
#                                        annotations,
#                                        by = "seurat_clusters")

load("figures/singleR_immgen_cell_labels.RData")
umap_df <- data.frame(UM1 = exp_agg@reductions$umap@cell.embeddings[,1],
                      UM2 = exp_agg@reductions$umap@cell.embeddings[,2],
                      cell_type = pred_main_labels)
umap_celltype_plot <- ggplot(umap_df, aes(x = UM1, y = UM2, color = cell_type)) + 
  geom_point(size = 0.5, alpha = 1) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  annotate("label", x = -4.5, y = 2.5, label = "Fibroblasts") + 
  annotate("label", x = 8.5, y = 10.5, label = "Endothelial Cells") + 
  annotate("label", x = 10.5, y = -8, label = "Epithelial Cells") +
  annotate("label",x = -2.8,y = -16.5,label = "B Cells") + 
  annotate("label", x = -1.5, y = -11.5, label = "T Cells") +
  annotate("label", x = 2, y = -5.9, label = "NKT Cells") + 
  annotate("label", x = -6.5, y = -10.5, label = "ILCs & NK Cells") +
  annotate("label", x = 2, y = -14, label = "NK Cells") + 
  annotate("label", x = 7.5, y = -11.5, label = "Monocytes & Macrophages")
ggsave(filename = "figures/umap_celltype_plot.pdf",
       plot = umap_celltype_plot,
       device = "pdf",
       width = 8,
       height = 7)
  
t_cell_df <- exp_agg@meta.data %>%
  mutate(time = as.numeric(week),
         cell_type = pred_main_labels) %>%
  group_by(time,group) %>%
  summarize(n_tcells = sum(cell_type == "T cells"))
b_cell_df <- exp_agg@meta.data %>%
  mutate(time = as.numeric(week),
         cell_type = pred_main_labels) %>%
  group_by(time,group) %>%
  summarize(n_bcells = sum(cell_type == "B cells"))
n_cell_df <- exp_agg@meta.data %>%
  mutate(time = as.numeric(week)) %>%
  group_by(time,group) %>%
  summarize(n_cells = n())
t_cell_df <- inner_join(t_cell_df,n_cell_df,by = c("time","group")) %>%
  mutate(prop_tcells = n_tcells/n_cells,
         se = 1.96*sqrt(prop_tcells*(1-prop_tcells)/n_cells),
         lb = prop_tcells - se,
         ub = prop_tcells + se)
b_cell_df <- inner_join(b_cell_df,n_cell_df,by = c("time","group")) %>%
  mutate(prop_bcells = n_bcells/n_cells,
         se = 1.96*sqrt(prop_bcells*(1-prop_bcells)/n_cells),
         lb = prop_bcells - se,
         ub = prop_bcells + se)

t_cell_prop_plot <- ggplot(data = t_cell_df, aes(x = time,
                                                 y = prop_tcells,
                                                 color = group,
                                                 group = group)) + 
  geom_point(size = 2) + 
  geom_line(size = 1.22) + 
  geom_errorbar(aes(ymin = lb, ymax = ub), size = 1) + 
  theme_classic() + 
  scale_x_continuous(breaks = c(10,10.5,17.5,23)) + 
  theme(legend.title = element_blank()) + 
  xlab("Weeks") + 
  ylab("Proportion of cells that are T cells (95% CI)") +
  theme(text = element_text(size = 15),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 315))
ggsave(t_cell_prop_plot,
       filename = "figures/t_cell_prop_plot.pdf",
       height = 6,
       width = 8)

b_cell_prop_plot <- ggplot(data = b_cell_df, aes(x = time,
                                                 y = prop_bcells,
                                                 color = group,
                                                 group = group)) + 
  geom_point(size = 2) + 
  geom_line(size = 1.22) + 
  geom_errorbar(aes(ymin = lb, ymax = ub), size = 1) + 
  theme_classic() + 
  scale_x_continuous(breaks = c(10,10.5,17.5,23)) + 
  theme(legend.title = element_blank()) + 
  xlab("Weeks") + 
  ylab("Proportion of cells that are B cells (95% CI)") +
  theme(text = element_text(size = 15),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 315))
ggsave(b_cell_prop_plot,
       filename = "figures/b_cell_prop_plot.pdf",
       height = 6,
       width = 8)

load("figures/singleR_immgen_cell_labels_fine.RData")
exp_agg@meta.data$cell_type_main <- pred_main_labels
exp_agg@meta.data$cell_type_fine <- pred_fine_labels
save(exp_agg,
     file = "data/Murine_BBN/integrated_processed_annotated_seurat.RData")
load("data/Murine_BBN/integrated_processed_annotated_seurat.RData")

AR <- exp_agg@assays$RNA@counts["Ar",]
AR_bin <- ifelse(AR == 0, "AR-","AR+")
exp_agg@meta.data$AR <- AR
exp_agg@meta.data$AR_bin <- AR_bin

umap_df$AR <- AR
umap_df$AR_bin <- AR_bin

ggplot(data = umap_df, aes(x = UM1, y = UM2, color = AR_bin)) + 
  geom_point()

exp_agg_tcells <- subset(exp_agg,
                         subset = cell_type_main == "T cells")
save(exp_agg_tcells,
     file = "data/Murine_BBN/exp_agg_tcells.RData")
exp_agg_imm <- subset(exp_agg,
                      subset = (cell_type_main != "Fibroblasts" &
                                cell_type_main != "Endothelial cells" &
                                cell_type_main != "Epithelial cells" &
                                cell_type_main != "Stromal cells" &
                                seurat_clusters != 7 &
                                seurat_clusters != 3 &
                                seurat_clusters != 14 &
                                seurat_clusters != 4))
save(exp_agg_imm,
     file = "data/Murine_BBN/exp_agg_immune.RData")

# DE analysis: control vs BBN
# Idents(exp_agg) <- "group"
# DE_BBN <- FindMarkers(exp_agg,
#                       ident.1 = "BBN",
#                       ident.2 = "control",
#                       logfc.threshold = 0,
#                       min.cells.group = 1,
#                       min.cells.feature = 1,
#                       min.pct = 0)
# write.csv(DE_BBN, file = "figures/DE_BBN/DE_BBN.csv")
# DE_BBN <- read.csv("figures/DE_BBN/DE_BBN.csv") %>%
#   mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
#          sig = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, TRUE, FALSE))
# top_DE_BBN_up <- DE_BBN %>%
#   filter(avg_log2FC > 0) %>%
#   top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
#   pull(X)
# top_DE_BBN_down <- DE_BBN %>%
#   filter(avg_log2FC < 0) %>%
#   top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
#   pull(X)
# label_genes_BBN <- c(top_DE_BBN_down,top_DE_BBN_up)
# v_BBN <- ggplot(data = DE_BBN,
#                  aes(x = -(avg_log2FC),
#                      y = -log10(p_val_adj + .Machine$double.xmin),
#                      color = sig)) + 
#   geom_point(alpha = 0.6) + 
#   scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) + 
#   geom_label_repel(aes(label = ifelse(X %in% label_genes_BBN,as.character(X),""),
#                        color = pol),
#                    box.padding   = 0.50, 
#                    point.padding = 0.30,
#                    label.size = 0.20) + 
#   theme_classic() + 
#   theme(legend.position = "none") + 
#   scale_y_continuous(expand = c(0.01,0.01)) +
#   xlab("Average log-fold change") + 
#   ylab("-log10(Adjusted p-value)")
# ggsave(v_BBN,
#        filename = "figures/BBN_volcano.pdf",
#        height = 6,
#        width = 8)


# packages
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(tidyverse)

# load counts & metadata
load('data/Murine_BBN/counts.RData')
# sce <- SingleCellExperiment(assays = list(counts = counts))
# 
# lognorm transform
# sce <- logNormCounts(sce)
# 
# load various reference databases
# ref_imm <- ImmGenData()
#
# cell type prediction
# pred_imm <- SingleR(test = sce, ref = ref_imm, labels = ref_imm$label.main)
# pred_imm_fine <- SingleR(test = sce, ref = ref_imm, labels = ref_imm$label.fine)
# 
# pred_list <- list(pred_imm,
#                   pred_imm_fine)
# 
# save(pred_list,
#      file = "data/Murine_BBN/singleR_predictions.RData")

load("data/Murine_BBN/singleR_predictions.RData")

pred_main <- pred_list[[1]]
pred_main_labels <- pred_main$labels
save(pred_main_labels,file = "figures/singleR_immgen_cell_labels.RData")

pred_fine <- pred_list[[2]]
pred_fine_labels <- pred_fine$labels
save(pred_fine_labels,file = "figures/singleR_immgen_cell_labels_fine.RData")

# prediction diagnostics - ImmGen Main 
pred_plot_imm <- plotScoreHeatmap(pred_list[[1]],
                                  annotation_col=as.data.frame(metadata[,"seurat_clusters",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_imm,
       filename = "figures/celltype_pred_heatmap_immgen.pdf",
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_imm <- table(pred_list[[1]]$labels, metadata[,"seurat_clusters"])
write.csv(pred_table_imm,
          file = "figures/celltype_pred_table_immgen.csv",
          quote = FALSE,
          row.names = TRUE)

# prediction diagnostics - ImmGen Fine 
pred_plot_imm_fine <- plotScoreHeatmap(pred_list[[2]],
                                       annotation_col=as.data.frame(metadata[,"seurat_clusters",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_imm_fine,
       filename = "figures/celltype_pred_heatmap_immgen_fine.pdf",
       device = "pdf",
       height  = 10,
       width  = 9)

pred_table_imm_fine <- table(pred_list[[2]]$labels, metadata[,"seurat_clusters"])
write.csv(pred_table_imm_fine,
          file = "figures/celltype_pred_table_immgen_fine.csv",
          quote = FALSE,
          row.names = TRUE)


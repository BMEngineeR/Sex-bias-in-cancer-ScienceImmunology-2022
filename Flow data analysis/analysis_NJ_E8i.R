# packages
library(data.table)
library(tidyverse)
library(ggpmisc)
library(ggrepel)
library(patchwork)
library(gmodels)
library(ggpubr)
library(rstatix)
library(hrbrthemes)
library(lme4)
library(lmerTest)
library(UCell)
library(FlowSOM)
library(premessa)
library(pheatmap)
# library(Rclusterpp)

# set seed
set.seed(1995)

data_dir <- "data/E8i_NJ/"
fig_dir <- "figures/E8i_NJ/"

# colors 
group_colors <- c("#7e0f5a","#006006")

# read gated data from NJ
# flow_df <- fread(paste0(data_dir,"18_files_concat-246727474865951_37.csv")) %>%
#   separate(col = "OmiqFileIndex",into = c("sex","ID"),sep = " ",remove = F) %>%
#   mutate(group = as.factor(substr(ID,1,2)))
# save(flow_df,
#      file = paste0(data_dir,"flow_df.RData"))
load(paste0(data_dir,"flow_df.RData"))

# use all markers for NJ's data
markers <- c("KI67",
             "CD27",
             "GITR",
             "LAG3",
             "CD62L",
             "ICOS",
             "CD95",
             "KLRG1",
             "VISTA",
             "TIGIT",
             "TIM3",
             "CD38",
             "TBET",
             "PD1",
             "EOMES",
             "TOX",
             "CTLA4",
             "CD69",
             "TCF1",
             "SLAMF6",
             "BCL2",
             "GZMB")
flow_mat <- flow_df %>%
  select(one_of(markers)) %>%
  as.matrix()


# umap_df <- uwot::umap(flow_mat)
# save(umap_df, file = paste0(data_dir,"umap_df.RData"))
load(paste0(data_dir,"umap_df.RData"))
umap_df <- as.data.frame(umap_df)
colnames(umap_df) <- c("UM1","UM2")

flow_df <- cbind(flow_df,umap_df) 
flow_df$day <- "D16"
metadata <- data.frame(cell = 1:nrow(flow_df),
                       stage = flow_df$day)

# PCA_gated <- princomp(flow_mat)
# save(PCA_gated,
#      file = paste0(data_dir,"PCA_gated.RData"))
load(paste0(data_dir,"PCA_gated.RData"))
PCA_df <- data.frame(PC1 = PCA_gated$scores[,1],
                     PC2 = PCA_gated$scores[,2],
                     sample = flow_df$OmiqFileIndex,
                     sex = flow_df$sex,
                     group = flow_df$group)

PCA_means <- PCA_df %>%
  group_by(sample,group) %>%
  summarize(PC1 = mean(PC1),
            PC2 = mean(PC2)) %>%
  ungroup()
PCA_sample_plot <- ggplot(PCA_means, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_text_repel(aes(label = sample)) +
  scale_color_manual(values = group_colors) +
  ggtitle("Location of samples on PC1-PC2 space") +
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(PCA_sample_plot,
       filename = paste0(fig_dir,"PCA_sample_plot.pdf"),
       height = 6,
       width = 9)

PCA_sample_plot_noID <- ggplot(PCA_means, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_manual(values = group_colors) +
  ggtitle("Location of samples on PC1-PC2 space") +
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(PCA_sample_plot_noID,
       filename = paste0(fig_dir,"PCA_sample_plot_noID.pdf"),
       height = 6,
       width = 9)

#######################################################################################################
# flowSOM workflow
f_frame <- premessa::as_flowFrame(flow_mat)
f_som <- FlowSOM(f_frame,
                 nClus = 13,
                 scale = FALSE,
                 xdim = 7,
                 ydim = 7)
save(f_som,
     file = paste0(data_dir,"f_som.RData"))
load(paste0(data_dir,"f_som.RData"))

# Plot results
PlotStars(f_som,
          backgroundValues = f_som$metaclustering)

PlotFlowSOM(f_som,
            backgroundValues = f_som$metaclustering,
            view = "MST")

# Get metaclustering per cell
f_som_clustering <- GetMetaclusters(f_som)
f_som_clustering_hires <- GetClusters(f_som)

# flow_som heatmap
flow_df$f_clust <- f_som_clustering
heatmap_df <- flow_df %>%
  select(-c(OmiqFileIndex,day,ID,sex,UM1,UM2,group))

markers_all <- colnames(heatmap_df)[1:ncol(heatmap_df)-1]
# markers_extra <- setdiff(markers_all,markers_ord)
# markers_ord_all <- c(markers_ord,markers_extra)

heatmap_df_long <- heatmap_df %>%
  pivot_longer(cols = -f_clust,
               names_to = "marker",
               values_to = "expression") %>%
  group_by(marker) %>%
  mutate(expression = (expression - median(expression))/(max(expression) - min(expression))) %>%
  ungroup()

# heatmap_df_sum <- heatmap_df_long %>%
#   group_by(marker,f_clust) %>%
#   summarize(expression = median(expression)) %>%
#   ungroup() %>%
#   mutate(marker = factor(marker,levels = markers_all))
# 
# heatmap_f_clust_all <- ggplot(data = heatmap_df_sum, aes(x = f_clust, y = marker, fill = expression)) +
#   geom_tile() +
#   # scale_fill_viridis_c(option = "A") +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   theme_classic() +
#   scale_x_discrete(expand = c(0,0))+
#   scale_y_discrete(expand = c(0,0)) + 
#   theme(legend.title = element_blank()) + 
#   ylab(NULL)
# ggsave(heatmap_cluster_all,
#        filename = paste0(fig_dir,"heatmap_f_clust_allmarkers.pdf"),
#        height = 5,
#        width = 7)

heatmap_df_long2 <- heatmap_df %>%
  pivot_longer(cols = -f_clust,
               names_to = "marker",
               values_to = "expression")

heatmap_df_sum2 <- heatmap_df_long2 %>%
  group_by(marker,f_clust) %>%
  summarize(expression = median(expression)) %>%
  ungroup() %>%
  pivot_wider(names_from = "f_clust",
              values_from = "expression")
mat_names <- heatmap_df_sum2$marker
heatmap_mat <- heatmap_df_sum2 %>%
  select(-marker) %>%
  as.matrix()
rownames(heatmap_mat) <- mat_names


heatmap_cluster_all2 <- pheatmap(heatmap_mat, 
                                 color = colorRampPalette(c("blue", "white","red"))(100),
                                 scale = "row", border_color = NA)
ggsave(heatmap_cluster_all2,
       filename = paste0(fig_dir,"heatmap_f_clust_allmarkers2.pdf"),
       height = 5,
       width = 7)

########################################################################################################
# Umap cluster labels
########################################################################################################

# CHOOSE ONE (cytotree or flowsom for the rest of the analysis)
umap_df$f_clust <- f_som_clustering
umap_df$cluster <- f_som_clustering

umap_df$sex <- flow_df$sex
umap_df$group <- flow_df$group
umap_df$day <- flow_df$day

umap_f_cluster_plot <- ggplot(data = umap_df, aes(x = UM1, y = UM2, color = f_clust)) + 
  geom_point(size = 0.1) + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  facet_grid(~ group)
ggsave(umap_f_cluster_plot,
       filename = paste0(fig_dir,"umap_f_cluster_plot.pdf"),
       height = 12,
       width = 12)

#########################################################################################################
# Umap marker expression plots
#########################################################################################################

umap_marker_CD69 <- umap_df %>%
  mutate(CD69 = flow_df$CD69) %>%
  mutate(CD69 = ifelse(CD69 < 0, 0, CD69)) %>%
  filter(CD69 < quantile(CD69,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (CD69))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("CD69") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD69,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","CD69.pdf"),
       height = 6,width = 6)

umap_marker_CD69_hex <- umap_df %>%
  mutate(CD69 = flow_df$CD69) %>%
  mutate(CD69 = ifelse(CD69 < 0, 0, CD69)) %>%
  filter(CD69 < quantile(CD69,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (CD69))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("CD69") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD69_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","CD69.pdf"),
       height = 6,width = 6)

umap_marker_CD27 <- umap_df %>%
  mutate(CD27 = flow_df$CD27) %>%
  mutate(CD27 = ifelse(CD27 < 0, 0, CD27)) %>%
  filter(CD27 < quantile(CD27,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (CD27))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("CD27") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD27,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","CD27.pdf"),
       height = 6,width = 6)

umap_marker_CD27_hex <- umap_df %>%
  mutate(CD27 = flow_df$CD27) %>%
  mutate(CD27 = ifelse(CD27 < 0, 0, CD27)) %>%
  filter(CD27 < quantile(CD27,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (CD27))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("CD27") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD27_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","CD27.pdf"),
       height = 6,width = 6)

umap_marker_LAG3 <- umap_df %>%
  mutate(LAG3 = flow_df$LAG3) %>%
  mutate(LAG3 = ifelse(LAG3 < 0, 0, LAG3)) %>%
  filter(LAG3 < quantile(LAG3,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (LAG3))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("LAG3") +
  theme(legend.title = element_blank())
ggsave(umap_marker_LAG3,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","LAG3.pdf"),
       height = 6,width = 6)

umap_marker_LAG3_hex <- umap_df %>%
  mutate(LAG3 = flow_df$LAG3) %>%
  mutate(LAG3 = ifelse(LAG3 < 0, 0, LAG3)) %>%
  filter(LAG3 < quantile(LAG3,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (LAG3))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("LAG3") +
  theme(legend.title = element_blank())
ggsave(umap_marker_LAG3_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","LAG3.pdf"),
       height = 6,width = 6)

umap_marker_CD95 <- umap_df %>%
  mutate(CD95 = flow_df$CD95) %>%
  mutate(CD95 = ifelse(CD95 < 0, 0, CD95)) %>%
  filter(CD95 < quantile(CD95,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (CD95))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("CD95") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD95,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","CD95.pdf"),
       height = 6,width = 6)

umap_marker_CD95_hex <- umap_df %>%
  mutate(CD95 = flow_df$CD95) %>%
  mutate(CD95 = ifelse(CD95 < 0, 0, CD95)) %>%
  filter(CD95 < quantile(CD95,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (CD95))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("CD95") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD95_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","CD95.pdf"),
       height = 6,width = 6)

umap_marker_KLRG1 <- umap_df %>%
  mutate(KLRG1 = flow_df$KLRG1) %>%
  mutate(KLRG1 = ifelse(KLRG1 < 0, 0, KLRG1)) %>%
  filter(KLRG1 < quantile(KLRG1,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (KLRG1))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("KLRG1") +
  theme(legend.title = element_blank())
ggsave(umap_marker_KLRG1,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","KLRG1.pdf"),
       height = 6,width = 6)

umap_marker_KLRG1_hex <- umap_df %>%
  mutate(KLRG1 = flow_df$KLRG1) %>%
  mutate(KLRG1 = ifelse(KLRG1 < 0, 0, KLRG1)) %>%
  filter(KLRG1 < quantile(KLRG1,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (KLRG1))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("KLRG1") +
  theme(legend.title = element_blank())
ggsave(umap_marker_KLRG1_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","KLRG1.pdf"),
       height = 6,width = 6)

umap_marker_VISTA <- umap_df %>%
  mutate(VISTA = flow_df$VISTA) %>%
  mutate(VISTA = ifelse(VISTA < 0, 0, VISTA)) %>%
  filter(VISTA < quantile(VISTA,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (VISTA))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("VISTA") +
  theme(legend.title = element_blank())
ggsave(umap_marker_VISTA,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","VISTA.pdf"),
       height = 6,width = 6)

umap_marker_VISTA_hex <- umap_df %>%
  mutate(VISTA = flow_df$VISTA) %>%
  mutate(VISTA = ifelse(VISTA < 0, 0, VISTA)) %>%
  filter(VISTA < quantile(VISTA,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (VISTA))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("VISTA") +
  theme(legend.title = element_blank())
ggsave(umap_marker_VISTA_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","VISTA.pdf"),
       height = 6,width = 6)


umap_marker_TIGIT <- umap_df %>%
  mutate(TIGIT = flow_df$TIGIT) %>%
  mutate(TIGIT = ifelse(TIGIT < 0, 0, TIGIT)) %>%
  filter(TIGIT < quantile(TIGIT,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (TIGIT))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("TIGIT") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TIGIT,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","TIGIT.pdf"),
       height = 6,width = 6)

umap_marker_TIGIT_hex <- umap_df %>%
  mutate(TIGIT = flow_df$TIGIT) %>%
  mutate(TIGIT = ifelse(TIGIT < 0, 0, TIGIT)) %>%
  filter(TIGIT < quantile(TIGIT,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (TIGIT))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("TIGIT") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TIGIT_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","TIGIT.pdf"),
       height = 6,width = 6)

umap_marker_TIM3 <- umap_df %>%
  mutate(TIM3 = flow_df$TIM3) %>%
  mutate(TIM3 = ifelse(TIM3 < 0, 0, TIM3)) %>%
  filter(TIM3 < quantile(TIM3,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (TIM3))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("TIM3") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TIM3,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","TIM3.pdf"),
       height = 6,width = 6)

umap_marker_TIM3_hex <- umap_df %>%
  mutate(TIM3 = flow_df$TIM3) %>%
  mutate(TIM3 = ifelse(TIM3 < 0, 0, TIM3)) %>%
  filter(TIM3 < quantile(TIM3,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (TIM3))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("TIM3") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TIM3_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","TIM3.pdf"),
       height = 6,width = 6)

umap_marker_TBET <- umap_df %>%
  mutate(TBET = flow_df$TBET) %>%
  mutate(TBET = ifelse(TBET < 0, 0, TBET)) %>%
  filter(TBET < quantile(TBET,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (TBET))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("TBET") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TBET,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","TBET.pdf"),
       height = 6,width = 6)

umap_marker_TBET_hex <- umap_df %>%
  mutate(TBET = flow_df$TBET) %>%
  mutate(TBET = ifelse(TBET < 0, 0, TBET)) %>%
  filter(TBET < quantile(TBET,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (TBET))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("TBET") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TBET_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","TBET.pdf"),
       height = 6,width = 6)


umap_marker_PD1 <- umap_df %>%
  mutate(PD1 = flow_df$PD1) %>%
  mutate(PD1 = ifelse(PD1 < 0, 0, PD1)) %>%
  filter(PD1 < quantile(PD1,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (PD1))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("PD1") +
  theme(legend.title = element_blank())
ggsave(umap_marker_PD1,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","PD1.pdf"),
       height = 6,width = 6)

umap_marker_PD1_hex <- umap_df %>%
  mutate(PD1 = flow_df$PD1) %>%
  mutate(PD1 = ifelse(PD1 < 0, 0, PD1)) %>%
  filter(PD1 < quantile(PD1,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (PD1))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("PD1") +
  theme(legend.title = element_blank())
ggsave(umap_marker_PD1_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","PD1.pdf"),
       height = 6,width = 6)


umap_marker_EOMES <- umap_df %>%
  mutate(EOMES = flow_df$EOMES) %>%
  mutate(EOMES = ifelse(EOMES < 0, 0, EOMES)) %>%
  filter(EOMES < quantile(EOMES,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (EOMES))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("EOMES") +
  theme(legend.title = element_blank())
ggsave(umap_marker_EOMES,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","EOMES.pdf"),
       height = 6,width = 6)

umap_marker_EOMES_hex <- umap_df %>%
  mutate(EOMES = flow_df$EOMES) %>%
  mutate(EOMES = ifelse(EOMES < 0, 0, EOMES)) %>%
  filter(EOMES < quantile(EOMES,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (EOMES))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("EOMES") +
  theme(legend.title = element_blank())
ggsave(umap_marker_EOMES_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","EOMES.pdf"),
       height = 6,width = 6)

umap_marker_TOX <- umap_df %>%
  mutate(TOX = flow_df$TOX) %>%
  mutate(TOX = ifelse(TOX < 0, 0, TOX)) %>%
  filter(TOX < quantile(TOX,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (TOX))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("TOX") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TOX,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","TOX.pdf"),
       height = 6,width = 6)

umap_marker_TOX_hex <- umap_df %>%
  mutate(TOX = flow_df$TOX) %>%
  mutate(TOX = ifelse(TOX < 0, 0, TOX)) %>%
  filter(TOX < quantile(TOX,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (TOX))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("TOX") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TOX_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","TOX.pdf"),
       height = 6,width = 6)

umap_marker_CTLA4 <- umap_df %>%
  mutate(CTLA4 = flow_df$CTLA4) %>%
  mutate(CTLA4 = ifelse(CTLA4 < 0, 0, CTLA4)) %>%
  filter(CTLA4 < quantile(CTLA4,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (CTLA4))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("CTLA4") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CTLA4,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","CTLA4.pdf"),
       height = 6,width = 6)

umap_marker_CTLA4_hex <- umap_df %>%
  mutate(CTLA4 = flow_df$CTLA4) %>%
  mutate(CTLA4 = ifelse(CTLA4 < 0, 0, CTLA4)) %>%
  filter(CTLA4 < quantile(CTLA4,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (CTLA4))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("CTLA4") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CTLA4_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","CTLA4.pdf"),
       height = 6,width = 6)

umap_marker_BCL2 <- umap_df %>%
  mutate(BCL2 = flow_df$BCL2) %>%
  mutate(BCL2 = ifelse(BCL2 < 0, 0, BCL2)) %>%
  filter(BCL2 < quantile(BCL2,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (BCL2))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("BCL2") +
  theme(legend.title = element_blank())
ggsave(umap_marker_BCL2,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","BCL2.pdf"),
       height = 6,width = 6)

umap_marker_BCL2_hex <- umap_df %>%
  mutate(BCL2 = flow_df$BCL2) %>%
  mutate(BCL2 = ifelse(BCL2 < 0, 0, BCL2)) %>%
  filter(BCL2 < quantile(BCL2,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (BCL2))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("BCL2") +
  theme(legend.title = element_blank())
ggsave(umap_marker_BCL2_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","BCL2.pdf"),
       height = 6,width = 6)

umap_marker_GZMB <- umap_df %>%
  mutate(GZMB = flow_df$GZMB) %>%
  mutate(GZMB = ifelse(GZMB < 0, 0, GZMB)) %>%
  filter(GZMB < quantile(GZMB,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (GZMB))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("GZMB") +
  theme(legend.title = element_blank())
ggsave(umap_marker_GZMB,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","GZMB.pdf"),
       height = 6,width = 6)

umap_marker_GZMB_hex <- umap_df %>%
  mutate(GZMB = flow_df$GZMB) %>%
  mutate(GZMB = ifelse(GZMB < 0, 0, GZMB)) %>%
  filter(GZMB < quantile(GZMB,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (GZMB))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("GZMB") +
  theme(legend.title = element_blank())
ggsave(umap_marker_GZMB_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","GZMB.pdf"),
       height = 6,width = 6)

umap_marker_GITR <- umap_df %>%
  mutate(GITR = flow_df$GITR) %>%
  mutate(GITR = ifelse(GITR < 0, 0, GITR)) %>%
  filter(GITR < quantile(GITR,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (GITR))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("GITR") +
  theme(legend.title = element_blank())
ggsave(umap_marker_GITR,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","GITR.pdf"),
       height = 6,width = 6)

umap_marker_GITR_hex <- umap_df %>%
  mutate(GITR = flow_df$GITR) %>%
  mutate(GITR = ifelse(GITR < 0, 0, GITR)) %>%
  filter(GITR < quantile(GITR,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (GITR))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("GITR") +
  theme(legend.title = element_blank())
ggsave(umap_marker_GITR_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","GITR.pdf"),
       height = 6,width = 6)

umap_marker_CD38 <- umap_df %>%
  mutate(CD38 = flow_df$CD38) %>%
  mutate(CD38 = ifelse(CD38 < 0, 0, CD38)) %>%
  filter(CD38 < quantile(CD38,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (CD38))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("CD38") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD38,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","CD38.pdf"),
       height = 6,width = 6)

umap_marker_CD38_hex <- umap_df %>%
  mutate(CD38 = flow_df$CD38) %>%
  mutate(CD38 = ifelse(CD38 < 0, 0, CD38)) %>%
  filter(CD38 < quantile(CD38,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (CD38))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("CD38") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD38_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","CD38.pdf"),
       height = 6,width = 6)

umap_marker_ICOS <- umap_df %>%
  mutate(ICOS = flow_df$ICOS) %>%
  mutate(ICOS = ifelse(ICOS < 0, 0, ICOS)) %>%
  filter(ICOS < quantile(ICOS,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (ICOS))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("ICOS") +
  theme(legend.title = element_blank())
ggsave(umap_marker_ICOS,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","ICOS.pdf"),
       height = 6,width = 6)

umap_marker_ICOS_hex <- umap_df %>%
  mutate(ICOS = flow_df$ICOS) %>%
  mutate(ICOS = ifelse(ICOS < 0, 0, ICOS)) %>%
  filter(ICOS < quantile(ICOS,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (ICOS))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("ICOS") +
  theme(legend.title = element_blank())
ggsave(umap_marker_ICOS_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","ICOS.pdf"),
       height = 6,width = 6)


umap_marker_CD62L <- umap_df %>%
  mutate(CD62L = flow_df$CD62L) %>%
  mutate(CD62L = ifelse(CD62L < 0, 0, CD62L)) %>%
  filter(CD62L < quantile(CD62L,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (CD62L))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("CD62L") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD62L,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","CD62L.pdf"),
       height = 6,width = 6)

umap_marker_CD62L_hex <- umap_df %>%
  mutate(CD62L = flow_df$CD62L) %>%
  mutate(CD62L = ifelse(CD62L < 0, 0, CD62L)) %>%
  filter(CD62L < quantile(CD62L,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (CD62L))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("CD62L") +
  theme(legend.title = element_blank())
ggsave(umap_marker_CD62L_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","CD62L.pdf"),
       height = 6,width = 6)

umap_marker_TCF1 <- umap_df %>%
  mutate(TCF1 = flow_df$TCF1) %>%
  mutate(TCF1 = ifelse(TCF1 < 0, 0, TCF1)) %>%
  filter(TCF1 < quantile(TCF1,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (TCF1))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("TCF1") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TCF1,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","TCF1.pdf"),
       height = 6,width = 6)

umap_marker_TCF1_hex <- umap_df %>%
  mutate(TCF1 = flow_df$TCF1) %>%
  mutate(TCF1 = ifelse(TCF1 < 0, 0, TCF1)) %>%
  filter(TCF1 < quantile(TCF1,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (TCF1))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("TCF1") +
  theme(legend.title = element_blank())
ggsave(umap_marker_TCF1_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","TCF1.pdf"),
       height = 6,width = 6)

umap_marker_SLAMF6 <- umap_df %>%
  mutate(SLAMF6 = flow_df$SLAMF6) %>%
  mutate(SLAMF6 = ifelse(SLAMF6 < 0, 0, SLAMF6)) %>%
  filter(SLAMF6 < quantile(SLAMF6,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (SLAMF6))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("SLAMF6") +
  theme(legend.title = element_blank())
ggsave(umap_marker_SLAMF6,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","SLAMF6.pdf"),
       height = 6,width = 6)

umap_marker_SLAMF6_hex <- umap_df %>%
  mutate(SLAMF6 = flow_df$SLAMF6) %>%
  mutate(SLAMF6 = ifelse(SLAMF6 < 0, 0, SLAMF6)) %>%
  filter(SLAMF6 < quantile(SLAMF6,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (SLAMF6))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("SLAMF6") +
  theme(legend.title = element_blank())
ggsave(umap_marker_SLAMF6_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","SLAMF6.pdf"),
       height = 6,width = 6)

umap_marker_KI67 <- umap_df %>%
  mutate(KI67 = flow_df$KI67) %>%
  mutate(KI67 = ifelse(KI67 < 0, 0, KI67)) %>%
  filter(KI67 < quantile(KI67,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, color = (KI67))) +
  geom_point(size = 0.1) +
  theme_void() +
  scale_color_viridis_c(option = "A") +
  ggtitle("KI67") +
  theme(legend.title = element_blank())
ggsave(umap_marker_KI67,
       filename = paste0(fig_dir,"marker_umap_plots/dot/","KI67.pdf"),
       height = 6,width = 6)

umap_marker_KI67_hex <- umap_df %>%
  mutate(KI67 = flow_df$KI67) %>%
  mutate(KI67 = ifelse(KI67 < 0, 0, KI67)) %>%
  filter(KI67 < quantile(KI67,probs = c(0.95))) %>%
  ggplot(data = .,
         aes(x = UM1, y = UM2, z = (KI67))) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  scale_fill_viridis_c(option = "A") +
  theme_void() +
  ggtitle("KI67") +
  theme(legend.title = element_blank())
ggsave(umap_marker_KI67_hex,
       filename = paste0(fig_dir,"marker_umap_plots/hex/","KI67.pdf"),
       height = 6,width = 6)


##############################################################################################################
# Signature scoring
#############################################################################################################

markers <- list()
markers$activation <- c("CD69+",
                        "GZMB+",
                        "GITR+",
                        "CD38+",
                        "ICOS+",
                        "KI67+")
#   CD62L-
#   Tcf1+
#   Slamf6+
#   Bcl2+
#   Tim3-
#   Lag3-
#   PD1-
#   Ctla4-
#   TOX-
markers$pe <- c("CD62L-",
                "TCF1+",
                "SLAMF6+",
                "TIM3-",
                "BCL2+",
                "LAG3-",
                "PD1-",
                "CTLA4-",
                "TOX-")
markers$naive <- c("CD62L+",
                   "TCF1+")
scores <- ScoreSignatures_UCell(t(flow_mat), features=markers)
save(scores,
     file = paste0(data_dir,"scores_ucell_all.RData"))
load(paste0(data_dir,"scores_ucell_all.RData"))
colnames(scores) <- c("Activation","PE","Naive")
umap_df <- cbind(umap_df,scores)

# umap score plots
# each score at each timepoint + overall
# remove outliers from visualizations
act_plot <- umap_df %>%
  filter(Activation < quantile(Activation,probs = 0.99) & Activation > quantile(Activation,probs = 0.01)) %>%
  ggplot(data = ., aes(x = UM1, y = UM2, color = Activation)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_viridis_c(option = "A") +
  ggtitle("Activation Signatures")
ggsave(act_plot,
       filename = paste0(fig_dir,"umap_activation.pdf"),
       height = 6,
       width = 9)

act_plot_hex <- umap_df %>%
  filter(Activation < quantile(Activation,probs = 0.99) & Activation > quantile(Activation,probs = 0.01)) %>%
  ggplot(data = ., aes(x = UM1, y = UM2, z = Activation)) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  theme_classic() +
  scale_fill_viridis_c(option = "A") +
  ggtitle("Activation Signatures")
ggsave(act_plot_hex,
       filename = paste0(fig_dir,"umap_activation_hex.pdf"),
       height = 6,
       width = 9)

pe_plot <- umap_df %>%
  # filter(PE < quantile(PE,probs = 0.99) & PE > quantile(PE,probs = 0.01)) %>%
  ggplot(data = ., aes(x = UM1, y = UM2, color = PE)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_viridis_c(option = "A") +
  ggtitle("Progenitor Exhaustion Signatures")
ggsave(pe_plot,
       filename = paste0(fig_dir,"umap_pe.pdf"),
       height = 6,
       width = 9)

pe_plot_hex <- umap_df %>%
  # filter(PE < quantile(PE,probs = 0.99) & PE > quantile(PE,probs = 0.01)) %>%
  ggplot(data = ., aes(x = UM1, y = UM2, z = PE)) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  theme_classic() +
  scale_fill_viridis_c(option = "A") +
  ggtitle("Progenitor Exhaustion Signatures")
ggsave(pe_plot_hex,
       filename = paste0(fig_dir,"umap_pe_hex.pdf"),
       height = 6,
       width = 9)

naive_plot <- umap_df %>%
  #filter(Naive < quantile(Naive,probs = 0.99)) %>%
  ggplot(data = ., aes(x = UM1, y = UM2, color = Naive)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_viridis_c(option = "A") +
  ggtitle("Naive T-cell Signatures")
ggsave(naive_plot,
       filename = paste0(fig_dir,"umap_naive.pdf"),
       height = 6,
       width = 9)

naive_plot_hex <- umap_df %>%
  #filter(Naive < quantile(Naive,probs = 0.99)) %>%
  ggplot(data = ., aes(x = UM1, y = UM2, z = Naive)) +
  stat_summary_hex(bins = 80,
                   color = "black",
                   size = 0.1) +
  theme_classic() +
  scale_fill_viridis_c(option = "A") +
  ggtitle("Naive T-cell Signatures")
ggsave(naive_plot_hex,
       filename = paste0(fig_dir,"umap_naive_hex.pdf"),
       height = 6,
       width = 9)

score_plots <- act_plot + pe_plot + naive_plot
ggsave(score_plots,
       filename = paste0(fig_dir,"umap_scores.pdf"),
       height = 6,
       width = 18)

score_plots_hex <- act_plot_hex + pe_plot_hex + naive_plot_hex
ggsave(score_plots_hex,
       filename = paste0(fig_dir,"umap_scores_hex.pdf"),
       height = 6,
       width = 18)

# Visualization for PE scores
#####################################################################################
# Cluster signature boxplots
#####################################################################################
cluster_PE_boxplots <- ggplot(data = umap_df,aes(x = reorder(cluster,-PE, FUN = median), y = PE)) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("MST Cluster") + 
  ylab("PE Signature Score") + 
  geom_hline(yintercept = quantile(umap_df$PE,probs = c(0.75)),
             linetype = "dashed") + 
  ggtitle("PE Signatures by MST Cluster") + 
  labs(subtitle = "Dashed Line Separates Non-PE Clusters")
ggsave(cluster_PE_boxplots, 
       filename = paste0(fig_dir,"cluster_PE_boxplots.pdf"),
       height = 6, 
       width = 9)

cluster_act_boxplots <- ggplot(data = umap_df,aes(x = reorder(cluster,-Activation), y = Activation)) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("MST Cluster") + 
  ylab("Activation Signature Score") + 
  geom_hline(yintercept = quantile(umap_df$Activation,probs = c(0.75)),
             linetype = "dashed") + 
  # geom_hline(yintercept = 0.999,
  #            linetype = "dashed") + 
  ggtitle("Activation Signatures by MST Cluster") + 
  labs(subtitle = "Dashed Line Separates Non-Activated Clusters")
ggsave(cluster_act_boxplots, 
       filename = paste0(fig_dir,"cluster_act_boxplots.pdf"),
       height = 6, 
       width = 9)

cluster_naive_boxplots <- ggplot(data = umap_df,aes(x = reorder(cluster,-Naive), y = Naive)) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab("MST Cluster") + 
  ylab("Naive Signature Score") + 
  geom_hline(yintercept = quantile(umap_df$Naive,probs = c(0.75)) + 0.0001,
             linetype = "dashed") + 
  ggtitle("Naive Signatures by MST Cluster") + 
  labs(subtitle = "Dashed Line Separates Non-Naive Clusters")
ggsave(cluster_naive_boxplots, 
       filename = paste0(fig_dir,"cluster_naive_boxplots.pdf"),
       height = 6, 
       width = 9)

################################################################################
# PE count expected plots - all time points
################################################################################
n_data <- nrow(flow_df)
umap_df$group <- umap_df$sex

PE_clusters <- umap_df %>%
  group_by(cluster) %>%
  summarize(med_PE = median(PE)) %>%
  #filter(med_PE > quantile(umap_df$PE,probs = c(0.75))) %>%
  top_n(n = 1,wt = med_PE) %>%
  ungroup() %>%
  pull(cluster)

umap_df <- umap_df %>%
  mutate(PE_cluster = ifelse(cluster %in% PE_clusters,1,0))
umap_df <- umap_df %>%
  mutate(PE_cell = ifelse(PE_cluster == 1,"PE Cell","Non PE Cell"))
group_counts <- umap_df %>%
  group_by(group) %>%
  summarize(n_group = n())
PE_count_df <- umap_df %>%
  group_by(PE_cell) %>%
  summarize(n_cell = n())
PE_group_counts <- umap_df %>%
  group_by(group,PE_cell) %>%
  summarize(n = n())
PE_group_counts_df <- inner_join(PE_group_counts,group_counts,by = "group")
PE_group_counts_df <- inner_join(PE_group_counts_df,PE_count_df,by = "PE_cell")
PE_group_counts_df <- PE_group_counts_df %>%
  ungroup() %>%
  mutate(n_data = n_data) 
PE_group_counts_df$e_count <- 0
for(i in 1:nrow(PE_group_counts_df))
{
  PE_group_counts_df$e_count[i] <- round(as.double((as.double(PE_group_counts_df$n_group[i])*as.double(PE_group_counts_df$n_cell[i]))/as.double(PE_group_counts_df$n_data[i])),0)
}
PE_group_counts_df$exp_count <- PE_group_counts_df$e_count
PE_group_counts_df <- PE_group_counts_df %>%
  mutate(adj_count = n-exp_count)
PE_count_sub <- PE_group_counts_df %>%
  filter(PE_cell == "PE Cell")
group_d <- umap_df %>%
  pull(group) %>%
  as.factor()
PE_d <- umap_df %>%
  pull(PE_cell) %>%
  as.factor()
chi_fit_d <- chisq.test(x = group_d,y = PE_d)
chi_fit_d$p.value
adjusted_PE_count_plot_all <- ggplot(PE_count_sub,aes(x = group,y = adj_count,fill = group)) +
  geom_col() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  ylab("Adjusted PE Cell Counts") +
  xlab(NULL) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-2*max(abs(PE_count_sub$adj_count)),2*max(abs(PE_count_sub$adj_count)))) +
  ggtitle("Day 7, 10, 13 combined") +
  labs(subtitle = paste0("p-value = ",round(chi_fit_d$p.value,3))) +
  annotate(geom = "text",x = 1, y = PE_count_sub$adj_count[1] + 10,label = paste0(abs(round((PE_count_sub$adj_count[1]/PE_count_sub$exp_count[1])*100,0)),"% vs. expected")) +
  annotate(geom = "text",x = 2, y = PE_count_sub$adj_count[2] - 10,label = paste0(abs(round((PE_count_sub$adj_count[2]/PE_count_sub$exp_count[2])*100,0)),"% vs. expected")) 
ggsave(adjusted_PE_count_plot_all,
       filename = paste0(fig_dir,"adjusted_PE_count_plot_all.pdf"),
       height = 6, 
       width = 8)

adjusted_PE_prop_plot_all <- 
  PE_count_sub %>%
  mutate(adj_perc = 100*(adj_count/exp_count)) %>%
  ggplot(data = .,aes(x = group,y = adj_perc,fill = group)) +
  geom_col() +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  ylab("% Difference Observed vs. Expected PE Cell Counts") +
  xlab(NULL) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-50,50)) +
  ggtitle("Day 7, 10, 13 combined")
ggsave(adjusted_PE_prop_plot_all,
       filename = paste0(fig_dir,"adjusted_PE_percentage_plot_all.pdf"),
       height = 6,
       width = 8)

################################################################################
# cell cluster count plots
################################################################################
umap_df$sample <- flow_df$ID

naive_clusters <- umap_df %>%
  group_by(cluster) %>%
  summarize(mean_naive = mean(Naive)) %>%
  top_n(n = 1,wt = mean_naive) %>%
  ungroup() %>%
  pull(cluster)

sample_counts <- umap_df %>%
  group_by(sample) %>%
  summarize(n_samp = n())

naive_sample_counts <- umap_df %>%
  filter(cluster %in% naive_clusters) %>%
  group_by(sample,group) %>%
  summarize(n_naive = n())
naive_sample_counts <- inner_join(naive_sample_counts,sample_counts,by = "sample") %>%
  mutate(prop_naive = n_naive/n_samp)
prop_naive_KO <- naive_sample_counts %>%
  filter(group == "KO") %>%
  pull(prop_naive)
prop_naive_WT <- naive_sample_counts %>%
  filter(group == "WT") %>%
  pull(prop_naive)

naive_sample_plots <- ggplot(naive_sample_counts, 
                             aes(x = group, y = prop_naive, fill = group)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1)+ 
  geom_jitter(size = 2,width = 0.1) + 
  theme_classic() + 
  scale_color_manual(values = group_colors) + 
  ggtitle("Naive Cell Proportions") + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_fill_manual(values = group_colors) + 
  theme(legend.position = "none") + 
  xlab(NULL) + 
  ylab("Cell Proportion")+
  labs(subtitle = paste0("Wilcoxon p-value = ",round(wilcox.test(prop_naive_KO,prop_naive_WT)$p.value,4)))


PE_clusters <- umap_df %>%
  group_by(cluster) %>%
  summarize(max_PE = max(PE)) %>%
  top_n(n = 1,wt = max_PE) %>%
  ungroup() %>%
  pull(cluster)

PE_sample_counts <- umap_df %>%
  filter(cluster %in% PE_clusters) %>%
  group_by(sample,group) %>%
  summarize(n_PE = n())
PE_sample_counts <- inner_join(PE_sample_counts,sample_counts,by = "sample") %>%
  mutate(prop_PE = n_PE/n_samp)
prop_PE_KO <- PE_sample_counts %>%
  filter(group == "KO") %>%
  pull(prop_PE)
prop_PE_WT <- PE_sample_counts %>%
  filter(group == "WT") %>%
  pull(prop_PE)

PE_sample_plots <- ggplot(PE_sample_counts, 
                          aes(x = group, y = prop_PE, fill = group)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1)+ 
  geom_jitter(size = 2,width = 0.1) + 
  theme_classic() + 
  scale_color_manual(values = group_colors) + 
  ggtitle("PE Cell Proporions") + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_fill_manual(values = group_colors) + 
  theme(legend.position = "none")+ 
  xlab(NULL)+ 
  ylab(NULL)+
  labs(subtitle = paste0("Wilcoxon p-value = ",round(wilcox.test(prop_PE_WT,prop_PE_KO)$p.value,4)))

act_clusters <- umap_df %>%
  group_by(cluster) %>%
  summarize(med_act = median(Activation)) %>%
  # filter(med_act > quantile(umap_df$Activation,probs = c(0.75))) %>%
  top_n(n = 1,wt = med_act) %>%
  ungroup() %>%
  pull(cluster)

act_sample_counts <- umap_df %>%
  filter(cluster %in% act_clusters) %>%
  group_by(sample,group) %>%
  summarize(n_act = n())
act_sample_counts <- inner_join(act_sample_counts,sample_counts,by = "sample") %>%
  mutate(prop_act = n_act/n_samp)
prop_act_KO <- act_sample_counts %>%
  filter(group == "KO") %>%
  pull(prop_act)
prop_act_WT <- act_sample_counts %>%
  filter(group == "WT") %>%
  pull(prop_act)

act_sample_plots <- ggplot(act_sample_counts, 
                           aes(x = group, y = prop_act, fill = group)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1)+ 
  geom_jitter(size = 2,width = 0.1) + 
  theme_classic() + 
  scale_color_manual(values = group_colors) + 
  ggtitle("Activated Cell Proportions") + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_fill_manual(values = group_colors) + 
  xlab(NULL)+ 
  ylab(NULL) + 
  labs(subtitle = paste0("Wilcoxon p-value = ",round(wilcox.test(prop_act_KO,prop_act_WT)$p.value,4)))

cell_prop_plots <- naive_sample_plots + PE_sample_plots + act_sample_plots
ggsave(cell_prop_plots,
       filename = paste0(fig_dir,"cell_prop_plots.pdf"),
       height = 6,
       width = 16)


# load packages
library(ggrepel)
library(tidyverse)
library(Seurat)
library(patchwork)

set.seed(1995)

# load immune cells
load("data/Murine_BBN/exp_agg_immune.RData")

# PD1 (PDCD1)?
PD1 <- exp_agg_imm@assays$RNA@counts["Pdcd1",]
PD1_bin <- ifelse(PD1 == 0,"PD1-","PD1+")

# SLAMF6
SLAMF6 <- exp_agg_imm@assays$RNA@counts["Slamf6",]
SLAMF6_bin <- ifelse(SLAMF6 == 0, "SLAMF6-", "SLAMF6+")

# TIM3 (HAVCR2)?
TIM3 <- exp_agg_imm@assays$RNA@counts["Havcr2",]
TIM3_bin <- ifelse(TIM3 == 0, "TIM3-", "TIM3+")

# Define progenitor exhausted cells
PE_status <- ifelse(PD1_bin == "PD1-" & SLAMF6_bin == "SLAMF6+" & TIM3_bin == "TIM3-",
                    "Progenitor exhausted",
                    "Not progenitor exhausted")

exp_agg_imm$PE <- PE_status

# load t cells
load("data/Murine_BBN/exp_agg_tcells.RData")
DimPlot(exp_agg_tcells,reduction = "umap")

DimPlot(exp_agg_imm,reduction = "umap")
umap_df_imm <- data.frame(UM1 = exp_agg_imm@reductions$umap@cell.embeddings[,1],
                          UM2 = exp_agg_imm@reductions$umap@cell.embeddings[,2],
                          cell_type_main = exp_agg_imm@meta.data$cell_type_main,
                          AR = exp_agg_imm@meta.data$AR,
                          TOX = exp_agg_imm@assays$RNA@counts["Tox",],
                          TCF7 = exp_agg_imm@assays$RNA@counts["Tcf7",],
                          SLAMF6 = exp_agg_imm@assays$RNA@counts["Slamf6",]) %>%
  mutate(TOX_bin = ifelse(TOX == 0, "TOX-", "TOX+"),
         TCF7_bin = ifelse(TCF7 == 0, "TCF7-", "TCF7+"),
         SLAMF6_bin = ifelse(SLAMF6 == 0, "SLAMF6-", "SLAMF6+"))
exp_agg_imm@meta.data$TOX_bin <- umap_df_imm$TOX_bin
exp_agg_imm@meta.data$TCF7_bin <- umap_df_imm$TCF7_bin
exp_agg_imm@meta.data$SLAMF6_bin <- umap_df_imm$SLAMF6_bin
write.csv(umap_df_imm,
          file = "immune_viz/data/umap_df_imm.csv",
          quote = FALSE,
          row.names = FALSE)
umap_celltype_plot_imm <- umap_df_imm %>%
  filter(UM2 > -17) %>% 
  ggplot(., aes(x = UM1, y = UM2, color = cell_type_main)) + 
  geom_point(size = 0.5, alpha = 1) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  annotate("label",x = -2.8,y = -16,label = "B Cells") + 
  annotate("label", x = -1.5, y = -11, label = "T Cells") +
  annotate("label", x = 2, y = -6.5, label = "NKT Cells") + 
  annotate("label", x = -6.75, y = -9.8, label = "NK Cells") + 
  annotate("label", x = -5.2, y = -8.3, label = "ILCs") + 
  annotate("label", x = 6, y = -11, label = "Monocytes") + 
  annotate("label", x = 6.1, y = -9.2, label = "Macrophages") + 
  annotate("label", x = -3, y = -7, label = "DCs") + 
  annotate("label", x = 2.5, y = -13, label = "ILCs") 
ggsave(filename = "figures/umap_celltype_plot_immune.pdf",
       plot = umap_celltype_plot_imm,
       device = "pdf",
       width = 8,
       height = 7)

ggplot(data = umap_df_imm, aes(x = UM1, y = UM2, color = AR)) + 
  geom_point(size = 0.8)

umap_AR_imm <- ggplot(data = umap_df_imm, aes(x = UM1, y = UM2, z = AR)) + 
  stat_summary_hex(bins = 50) + 
  scale_fill_viridis_c(option = "A") + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  annotate("label",x = -2.8,y = -16,label = "B Cells") + 
  annotate("label", x = -1.5, y = -11, label = "T Cells") +
  annotate("label", x = 2, y = -6.5, label = "NKT Cells") + 
  annotate("label", x = -6.75, y = -9.8, label = "NK Cells") + 
  annotate("label", x = -5, y = -8.3, label = "ILCs") + 
  annotate("label", x = 6, y = -11, label = "Monocytes") + 
  annotate("label", x = 6.1, y = -9.2, label = "Macrophages") + 
  annotate("label", x = -3, y = -7, label = "DCs") + 
  annotate("label", x = 2.5, y = -13, label = "ILCs")
ggsave(filename = "figures/umap_AR_immune.pdf",
       plot = umap_AR_imm,
       device = "pdf",
       width = 8,
       height = 7)

# DE AR positive vs negative in immune cells
Idents(exp_agg_imm) <- "AR_bin"
# DE_AR_immune <- FindMarkers(exp_agg_imm,
#                             ident.1 = "AR-",
#                             ident.2 = "AR+",
#                             logfc.threshold = 0,
#                             min.cells.group = 1,
#                             min.cells.feature = 1,
#                             min.pct = 0)
# write.csv(DE_AR_immune, file = "figures/DE_AR/DE_AR_immune.csv")
DE_AR_immune <- read.csv("figures/DE_AR/DE_AR_immune.csv") %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
         sig = ifelse(p_val_adj < 0.05, TRUE, FALSE))
top_DE_AR_immune_up <- DE_AR_immune %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_DE_AR_immune_down <- DE_AR_immune %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
label_genes_DE_AR_immune <- c(top_DE_AR_immune_up,top_DE_AR_immune_down)
v_AR_immune <- ggplot(data = DE_AR_immune,
                 aes(x = (avg_log2FC),
                     y = -log10(p_val_adj + .Machine$double.xmin),
                     color = sig)) +
  geom_point(alpha = 0.6) +
  annotate("label", x = 1.8, y = 18, label = "Up in AR-")+
  annotate("label", x = -1.8, y = 18, label = "Up in AR+") + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) +
  geom_label_repel(aes(label = ifelse(X %in% label_genes_DE_AR_immune,as.character(X),""),
                       color = pol),
                   box.padding   = 0.90,
                   point.padding = 0.70,
                   label.size = 0.15) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-2.5,2.5)) + 
  xlab("Average log-fold change") +
  ylab("-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0) 
ggsave(v_AR_immune,
       filename = "figures/AR_immune_volcano.pdf",
       height = 6,
       width = 8)

# DE TOX positive vs negative in immune cells
Idents(exp_agg_imm) <- "TOX_bin"
# DE_TOX_immune <- FindMarkers(exp_agg_imm,
#                             ident.1 = "TOX-",
#                             ident.2 = "TOX+",
#                             logfc.threshold = 0,
#                             min.cells.group = 1,
#                             min.cells.feature = 1,
#                             min.pct = 0)
# write.csv(DE_TOX_immune, file = "figures/DE_TOX/DE_TOX_immune.csv")
DE_TOX_immune <- read.csv("figures/DE_TOX/DE_TOX_immune.csv") %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
         sig = ifelse(p_val_adj < 0.05, TRUE, FALSE))
top_DE_TOX_immune_up <- DE_TOX_immune %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_DE_TOX_immune_down <- DE_TOX_immune %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
label_genes_DE_TOX_immune <- c(top_DE_TOX_immune_up,top_DE_TOX_immune_down)
v_TOX_immune <- ggplot(data = DE_TOX_immune,
                      aes(x = (avg_log2FC),
                          y = -log10(p_val_adj + .Machine$double.xmin),
                          color = sig)) +
  geom_point(alpha = 0.6) +
  annotate("label", x = 1.8, y = 25, label = "Up in TOX-")+
  annotate("label", x = -1.8, y = 25, label = "Up in TOX+") + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) +
  geom_label_repel(aes(label = ifelse(X %in% label_genes_DE_TOX_immune,as.character(X),""),
                       color = pol),
                   box.padding   = 0.90,
                   point.padding = 0.70,
                   label.size = 0.15,
                   max.overlaps = Inf) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-2.5,2.5)) + 
  xlab("Average log-fold change") +
  ylab("-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0) 
ggsave(v_TOX_immune,
       filename = "figures/TOX_immune_volcano.pdf",
       height = 6,
       width = 8)


# DE TCF7 positive vs negative in immune cells
Idents(exp_agg_imm) <- "TCF7_bin"
# DE_TCF7_immune <- FindMarkers(exp_agg_imm,
#                              ident.1 = "TCF7-",
#                              ident.2 = "TCF7+",
#                              logfc.threshold = 0,
#                              min.cells.group = 1,
#                              min.cells.feature = 1,
#                              min.pct = 0)
# write.csv(DE_TCF7_immune, file = "figures/DE_TCF7/DE_TCF7_immune.csv")
DE_TCF7_immune <- read.csv("figures/DE_TCF7/DE_TCF7_immune.csv") %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
         sig = ifelse(p_val_adj < 0.05, TRUE, FALSE)) %>%
  filter(X != "Tcf7")
top_DE_TCF7_immune_up <- DE_TCF7_immune %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_DE_TCF7_immune_down <- DE_TCF7_immune %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
label_genes_DE_TCF7_immune <- c(top_DE_TCF7_immune_up,top_DE_TCF7_immune_down)
v_TCF7_immune <- ggplot(data = DE_TCF7_immune,
                       aes(x = (avg_log2FC),
                           y = -log10(p_val_adj + .Machine$double.xmin),
                           color = sig)) +
  geom_point(alpha = 0.6) +
  annotate("label", x = 1.8, y = 50, label = "Up in TCF7-")+
  annotate("label", x = -1.8, y = 50, label = "Up in TCF7+") + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) +
  geom_label_repel(aes(label = ifelse(X %in% label_genes_DE_TCF7_immune,as.character(X),""),
                       color = pol),
                   box.padding   = 0.90,
                   point.padding = 0.70,
                   label.size = 0.15,
                   max.overlaps = Inf) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-2.5,2.5)) + 
  xlab("Average log-fold change") +
  ylab("-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0) 
ggsave(v_TCF7_immune,
       filename = "figures/TCF7_immune_volcano.pdf",
       height = 6,
       width = 8)



# DE SLAMF6 positive vs negative in immune cells
Idents(exp_agg_imm) <- "SLAMF6_bin"
# DE_SLAMF6_immune <- FindMarkers(exp_agg_imm,
#                               ident.1 = "SLAMF6-",
#                               ident.2 = "SLAMF6+",
#                               logfc.threshold = 0,
#                               min.cells.group = 1,
#                               min.cells.feature = 1,
#                               min.pct = 0)
# write.csv(DE_SLAMF6_immune, file = "figures/DE_SLAMF6/DE_SLAMF6_immune.csv")
DE_SLAMF6_immune <- read.csv("figures/DE_SLAMF6/DE_SLAMF6_immune.csv") %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
         sig = ifelse(p_val_adj < 0.05, TRUE, FALSE)) %>%
  filter(X != "Slamf6")
top_DE_SLAMF6_immune_up <- DE_SLAMF6_immune %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_DE_SLAMF6_immune_down <- DE_SLAMF6_immune %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
label_genes_DE_SLAMF6_immune <- c(top_DE_SLAMF6_immune_up,top_DE_SLAMF6_immune_down)
v_SLAMF6_immune <- ggplot(data = DE_SLAMF6_immune,
                        aes(x = (avg_log2FC),
                            y = -log10(p_val_adj + .Machine$double.xmin),
                            color = sig)) +
  geom_point(alpha = 0.6) +
  annotate("label", x = 1.8, y = 15, label = "Up in SLAMF6-")+
  annotate("label", x = -1.8, y = 15, label = "Up in SLAMF6+") + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) +
  geom_label_repel(aes(label = ifelse(X %in% label_genes_DE_SLAMF6_immune,as.character(X),""),
                       color = pol),
                   box.padding   = 0.90,
                   point.padding = 0.70,
                   label.size = 0.15,
                   max.overlaps = Inf) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-2.5,2.5)) + 
  xlab("Average log-fold change") +
  ylab("-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0) 
ggsave(v_SLAMF6_immune,
       filename = "figures/SLAMF6_immune_volcano.pdf",
       height = 6,
       width = 8)


# DE AR positive vs negative in T cells
Idents(exp_agg_tcells) <- "AR_bin"
# DE_AR_tcells <- FindMarkers(exp_agg_tcells,
#                             ident.1 = "AR-",
#                             ident.2 = "AR+",
#                             logfc.threshold = 0,
#                             min.cells.group = 1,
#                             min.cells.feature = 1,
#                             min.pct = 0)
# write.csv(DE_AR_tcells, file = "figures/DE_AR/DE_AR_tcells.csv")
DE_AR_tcells <- read.csv("figures/DE_AR/DE_AR_tcells.csv") %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
         sig = ifelse(p_val_adj < 0.05, TRUE, FALSE))
top_DE_AR_tcells_up <- DE_AR_tcells %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_DE_AR_tcells_down <- DE_AR_tcells %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
label_genes_DE_AR_tcells <- c(top_DE_AR_tcells_up,top_DE_AR_tcells_down)
v_AR_tcells <- ggplot(data = DE_AR_tcells,
                      aes(x = -(avg_log2FC),
                          y = -log10(p_val_adj + .Machine$double.xmin),
                          color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) +
  geom_label_repel(aes(label = ifelse(X %in% label_genes_DE_AR_tcells,as.character(X),""),
                       color = pol),
                   box.padding   = 0.90,
                   point.padding = 0.70,
                   label.size = 0.15) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(-2.5,2.5)) + 
  xlab("Average log-fold change") +
  ylab("-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0)
ggsave(v_AR_tcells,
       filename = "figures/AR_tcells_volcano.pdf",
       height = 6,
       width = 8)



# Specific subset analysis
# From Dr. Sundi:
# A couple of directions that were indicated this week 
# were to compare Tcf7 and other differentially expressed 
# genes (among T cells) in specific ways:
#   
# 1) Progenitor exhausted vs others
# 2) Antigen experienced vs antigen naive
# 3) Over time

# Progenitor exhausted
# PD1-
# Slamf6+
# Tim3-

# PD1 (PDCD1)?
PD1 <- exp_agg_tcells@assays$RNA@counts["Pdcd1",]
PD1_bin <- ifelse(PD1 == 0,"PD1-","PD1+")

# SLAMF6
SLAMF6 <- exp_agg_tcells@assays$RNA@counts["Slamf6",]
SLAMF6_bin <- ifelse(SLAMF6 == 0, "SLAMF6-", "SLAMF6+")

# TIM3 (HAVCR2)?
TIM3 <- exp_agg_tcells@assays$RNA@counts["Havcr2",]
TIM3_bin <- ifelse(TIM3 == 0, "TIM3-", "TIM3+")

# Define progenitor exhausted cells
PE_status <- ifelse(PD1_bin == "PD1-" & SLAMF6_bin == "SLAMF6+" & TIM3_bin == "TIM3-",
                    "Progenitor exhausted",
                    "Not progenitor exhausted")

# Terminally exhausted
# PD1+
# Slamf6-
# Tim3+

# Define terminally exhausted cells
TE_status <- ifelse(PD1_bin == "PD1+" & SLAMF6_bin == "SLAMF6-" & TIM3_bin == "TIM3+",
                    "Terminally exhausted",
                    "Not terminally exhausted")

# Antigen experienced
# CD44+
# CD62L-

CD44 <- exp_agg_tcells@assays$RNA@counts["Cd44",]
CD44_bin <- ifelse(CD44 == 0,"CD44-","CD44+")

# Can't find 
# Sell
CD62L <- exp_agg_tcells@assays$RNA@counts["Sell",]
CD62L_bin <- ifelse(CD62L == 0,"CD62L-","CD62L+")

# Antigen naive
# CD44-
# CD62L+
AN_status <- ifelse(CD44_bin == "CD44-" & CD62L_bin == "CD62L+",
                    "Antigen naive",
                    "Not antigen naive")

# Antigen experienced
# CD44+
# CD62L-
AE_status <- ifelse(CD44_bin == "CD44+" & CD62L_bin == "CD62L-",
                    "Antigen experienced",
                    "Not antigen experienced")

umap_tcells <- data.frame(UM1 = exp_agg_tcells@reductions$umap@cell.embeddings[,1],
                          UM2 = exp_agg_tcells@reductions$umap@cell.embeddings[,2],
                          PE = PE_status,
                          TE = TE_status,
                          AN = AN_status,
                          AE = AE_status,
                          TCF7 = exp_agg_tcells@assays$RNA@counts["Tcf7",],
                          SLAMF6 = exp_agg_tcells@assays$RNA@counts["Slamf6",])

exp_agg_tcells$PE <- umap_tcells$PE

umap_tcells_PE <- ggplot(data = umap_tcells, aes(x = UM1, y = UM2, color = PE)) + 
  geom_point(size = 0.6) + 
  theme_classic() + 
  facet_wrap(~ PE)
ggsave(umap_tcells_PE,
       filename = "figures/umap_tcells_PE.pdf",
       height = 5,
       width = 6)

umap_tcells %>%
  group_by(PE) %>%
  summarize(c_val = cor(TCF7,SLAMF6),
            p_val = cor.test(TCF7,SLAMF6)$p.value)

umap_tcells_TE <- ggplot(data = umap_tcells, aes(x = UM1, y = UM2, color = TE)) + 
  geom_point(size = 0.6) + 
  theme_classic() + 
  facet_wrap(~ TE)
ggsave(umap_tcells_TE,
       filename = "figures/umap_tcells_TE.pdf",
       height = 5,
       width = 6)

umap_tcells_AN <- ggplot(data = umap_tcells, aes(x = UM1, y = UM2, color = AN)) + 
  geom_point(size = 0.6) + 
  theme_classic() + 
  facet_wrap(~ AN)
ggsave(umap_tcells_AN,
       filename = "figures/umap_tcells_AN.pdf",
       height = 5,
       width = 6)

umap_tcells_AE <- ggplot(data = umap_tcells, aes(x = UM1, y = UM2, color = AE)) + 
  geom_point(size = 0.6) + 
  theme_classic() + 
  facet_wrap(~ AE)
ggsave(umap_tcells_AE,
       filename = "figures/umap_tcells_AE.pdf",
       height = 5,
       width = 6)

tcf7_boxplot_PE <- ggplot(data = umap_tcells, aes(x = PE_status, y = TCF7)) + 
  geom_boxplot() + 
  theme_classic()
ggsave(tcf7_boxplot_PE,
       filename = "figures/tcf7_boxplot_PE.pdf",
       height = 5,
       width = 6)

tcf7_boxplot_TE <- ggplot(data = umap_tcells, aes(x = TE_status, y = TCF7)) + 
  geom_boxplot() + 
  theme_classic()
ggsave(tcf7_boxplot_TE,
       filename = "figures/tcf7_boxplot_TE.pdf",
       height = 5,
       width = 6)

tcf7_boxplot_AN <- ggplot(data = umap_tcells, aes(x = AN_status, y = TCF7)) + 
  geom_boxplot() + 
  theme_classic()
ggsave(tcf7_boxplot_AN,
       filename = "figures/tcf7_boxplot_AN.pdf",
       height = 5,
       width = 6)

tcf7_boxplot_AE <- ggplot(data = umap_tcells, aes(x = AE_status, y = TCF7)) + 
  geom_boxplot() + 
  theme_classic()
ggsave(tcf7_boxplot_AE,
       filename = "figures/tcf7_boxplot_AE.pdf",
       height = 5,
       width = 6)

tcf7_boxplot_ANAE <- umap_tcells %>%
  filter(AN == "Antigen naive" | AE == "Antigen experienced") %>%
  mutate(ANAE = ifelse(AN == "Antigen naive", "Antigen naive", "Antigen experienced")) %>%
  ggplot(data = ., aes(x = ANAE, y = TCF7)) + 
  geom_boxplot() + 
  theme_classic()
ggsave(tcf7_boxplot_ANAE,
       filename = "figures/tcf7_boxplot_ANAE.pdf",
       height = 5,
       width = 6)




# Sub cluster T & NKT cells 
table(exp_agg_imm$cell_type_main)
exp_agg_tcells2 <- subset(exp_agg_imm,
                          subset = cell_type_main %in% c("T cells","NKT"))

DefaultAssay(exp_agg_tcells2) <- "integrated"
DimPlot(exp_agg_tcells2,reduction = "umap")

# Run the standard workflow for visualization and clustering
exp_agg_tcells2 <- ScaleData(exp_agg_tcells2, verbose = FALSE)
exp_agg_tcells2 <- RunPCA(exp_agg_tcells2, npcs = 30, verbose = FALSE)

# umap and Clustering
exp_agg_tcells2 <- RunUMAP(exp_agg_tcells2, reduction = "pca", dims = 1:20)
exp_agg_tcells2 <- FindNeighbors(exp_agg_tcells2, reduction = "pca", dims = 1:20)
exp_agg_tcells2 <- FindClusters(exp_agg_tcells2, resolution = 0.2)

p1 <- Seurat::DimPlot(exp_agg_tcells2, reduction = "umap", label = TRUE) + 
  theme(legend.position = "none")
ggsave(p1, 
       filename = "figures/umap_t_nkt_cells.pdf",
       height = 6,
       width = 8)

g_vec <- c("Cd3e","Cd8a","Cd4")
p2 <- FeaturePlot(exp_agg_tcells2,features = g_vec)
p3 <- VlnPlot(exp_agg_tcells2,features = "Cd8a",pt.size = 0.1)

cd8_plot <- p1 + p3
ggsave(cd8_plot, filename = "figures/umap_cd8_violin.pdf",
       height = 4,
       width = 8)




DimPlot(exp_agg_tcells2,reduction = "umap")
# sub cluster CD8+ t cells
exp_agg_CD8 <- subset(exp_agg_tcells2, 
                      subset = seurat_clusters %in% c(0,2,5))

DefaultAssay(exp_agg_CD8) <- "integrated"
DimPlot(exp_agg_CD8,reduction = "umap",label = TRUE)

# Run the standard workflow for visualization and clustering
exp_agg_CD8 <- ScaleData(exp_agg_CD8, verbose = FALSE)
exp_agg_CD8 <- RunPCA(exp_agg_CD8, npcs = 30, verbose = FALSE)

# umap and Clustering
exp_agg_CD8 <- RunUMAP(exp_agg_CD8, reduction = "pca", dims = 1:20)
exp_agg_CD8 <- FindNeighbors(exp_agg_CD8, reduction = "pca", dims = 1:20)
exp_agg_CD8 <- FindClusters(exp_agg_CD8, resolution = 0.2)

save(exp_agg_CD8,
     file = "data/integrated_processed_CD8T.RData")

p1 <- Seurat::DimPlot(exp_agg_CD8, reduction = "umap", label = TRUE) + 
  theme(legend.position = "none")
ggsave(p1,
       filename = "figures/umap_CD8T.pdf",
       height = 4,
       width = 6)


umap_CD8 <- data.frame(UM1 = exp_agg_CD8@reductions$umap@cell.embeddings[,1],
                       UM2 = exp_agg_CD8@reductions$umap@cell.embeddings[,2],
                       PE = exp_agg_CD8$PE,
                       TCF7 = exp_agg_CD8@assays$RNA@counts["Tcf7",],
                       SLAMF6 = exp_agg_CD8@assays$RNA@counts["Slamf6",])

umap_CD8 %>%
  group_by(PE) %>%
  summarize(n = n(),
            c_val = cor(TCF7,SLAMF6),
            p_val = cor.test(TCF7,SLAMF6)$p.value)

umap_CD8 %>%
  summarize(n = n(),
            c_val = cor(TCF7,SLAMF6),
            p_val = cor.test(TCF7,SLAMF6)$p.value)


Idents(exp_agg_CD8) <- "seurat_clusters"
g2_vec <- c("SELL", 
            "TCF7", 
            "LEF1", 
            "IL7R", 
            "CCR7", 
            "BCL2", 
            "KLRG1", 
            "IL2", 
            "GZMA", 
            "GNLY", 
            "GZMB", 
            "GZMK", 
            "IFNG", 
            "NKG7", 
            "HAVCR2", 
            "PDCD1", 
            "LAG3", 
            "TIGIT", 
            "CTLA4", 
            "CD28", 
            "TNFRSF14", 
            "ICOS", 
            "TNFRSF9", 
            "EOMES", 
            "HOPX", 
            "TBX21", 
            "ZEB2", 
            "ZNF683", 
            "HIF1A", 
            "ID2", 
            "TOX")

g3_vec <- rownames(exp_agg_CD8@assays$RNA@counts)[which(toupper(rownames(exp_agg_CD8@assays$RNA@counts)) %in% g2_vec)]
DefaultAssay(exp_agg_CD8) <- "RNA"
DoHeatmap(exp_agg_CD8,features = g3_vec, slot = "counts")
VlnPlot(exp_agg_CD8,features = g3_vec)

g_mat <- t(as.matrix(exp_agg_CD8@assays$RNA@counts[g3_vec,]))
g_df <- as.data.frame(g_mat)
g_df$cluster <- exp_agg_CD8$seurat_clusters

g_df_long <- g_df %>%
  pivot_longer(cols = -"cluster",
               names_to = "gene",
               values_to = "count") %>%
  group_by(gene) %>%
  mutate(z_count = scale(count))

g_means <- g_df_long %>%
  group_by(cluster,gene) %>%
  summarize(mean_exp = mean(z_count)) %>%
  pivot_wider(names_from = "gene",
              values_from = "mean_exp")

g_mean_mat <- as.matrix(g_means[,-1])
rownames(g_mean_mat) <- 1:5
g_mean_mat <- scale(g_mean_mat)
g_heatmap <- pheatmap::pheatmap(g_mean_mat,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE)
ggsave(g_heatmap,
       filename = "figures/CD8T_heatmap.pdf",
       height = 5,
       width = 7)



# TODO:
# 1) Make heatmap above
# 2) Plot violin plot of marker genes with umap feature plot
# send to Tong

for(g in g3_vec)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,features = g,pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T/",g,".pdf"),
         height = 8,
         width = 6)
}

# Additional genes
# Apoptosis-related: Fas
# Activation/naïve state-related: Slamf6, Cd69, Cx3cr1
# Ar-related: Ar, Ifna5, Infb1, Igf1, Cxcl9, Il6, Fgf8, Fkbp5, Igfbp3, Klkb1, Klk4, Mmp2, Mme, Nkx3-1, Ccn5.

apop_genes <- c("Fas")
act_genes <- c("Slamf6","Cd69","Cx3cr1")
ar_genes <- c("Ar",
              #"Ifna5",
              #"Infb1",
              "Igf1",
              "Cxcl9",
              "Il6",
              #"Fgf8",
              "Fkbp5",
              "Igfbp3",
              #"Klkb1",
              #"Klk4",
              "Mmp2",
              "Mme",
              #"Nkx3-1",
              "Ccn5")

new_genes <- c("Cd3e"
               ,"Cd4"
               ,"Il2ra"
               ,"Cd44"
               ,"Prf1"
               ,"Entpd1"
               ,"Cxcl13"
               ,"Cd101"
               ,"Cd38"
               ,"Vsir"
               ,"Tnfrsf18"
               ,"Id3"
               ,"Irf4"
               ,"Batf"
               ,"Nfatc1"
               ,"Nfatc2"
               ,"Nr4a1"
               ,"Nr4a2"
               ,"Nr4a3")

re_genes <- read.csv("Genes_re-draw_BBN dataset.csv") %>%
  pull(genes)
for(g in re_genes)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,features = g,pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_genes_redraw/",g,".pdf"),
         height = 8,
         width = 6)
}

re_genes2 <- c("Cd8a",
               #"Gnly",
               "Klrb1c",
               "Ncam1")
for(g in re_genes2)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,features = g,pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_genes_redraw2/",g,".pdf"),
         height = 8,
         width = 6)
}

for(g in apop_genes)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,features = g,pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_apoptosis/",g,".pdf"),
         height = 8,
         width = 6)
}

for(g in apop_genes)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,features = g,pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_apoptosis/",g,".pdf"),
         height = 8,
         width = 6)
}

for(g in act_genes)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,features = g,pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_activation/",g,".pdf"),
         height = 8,
         width = 6)
}

for(g in ar_genes)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,
                      features = g,
                      pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,
                          features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_AR/",g,".pdf"),
         height = 8,
         width = 6)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_AR/",g,".jpg"),
         height = 8,
         width = 6)
}

for(g in new_genes)
{
  v_plot_g <- VlnPlot(exp_agg_CD8,
                      features = g,
                      pt.size = 0.5)
  f_plot_g <- FeaturePlot(exp_agg_CD8,
                          features = g)
  plot_g <- v_plot_g + f_plot_g + plot_layout(ncol = 1)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_newgenes/",g,".pdf"),
         height = 8,
         width = 6)
  ggsave(plot_g,
         filename = paste0("figures/CD8T_newgenes/",g,".jpg"),
         height = 8,
         width = 6)
}


cor(exp_agg_CD8@assays$integrated@scale.data["Tcf7",],exp_agg_CD8@assays$integrated@scale.data["Slamf6",])
cor.test(exp_agg_CD8@assays$integrated@scale.data["Tcf7",],
         exp_agg_CD8@assays$integrated@scale.data["Slamf6",],
         method = "spearman")
cor_df <- data.frame(Tcf7 = exp_agg_CD8@assays$integrated@scale.data["Tcf7",],
                     Slamf6 = exp_agg_CD8@assays$integrated@scale.data["Slamf6",])
cor_plot <- ggplot(cor_df, aes(x = Tcf7, y = Slamf6)) + 
  geom_point() + 
  theme_classic() 
ggsave(cor_plot,
       filename = "figures/cor_plot_CD8T_TCF7_SLAMF6.pdf",
       height = 6,
       width = 8)

exp_agg_CD8_PE1 <- subset(exp_agg_CD8,
                          subset = seurat_clusters %in% c(0,4))
exp_agg_CD8_PE2 <- subset(exp_agg_CD8,
                          subset = seurat_clusters %in% c(0))


cor(exp_agg_CD8_PE1@assays$integrated@scale.data["Tcf7",],exp_agg_CD8_PE1@assays$integrated@scale.data["Slamf6",])
cor.test(exp_agg_CD8_PE1@assays$integrated@scale.data["Tcf7",],
         exp_agg_CD8_PE1@assays$integrated@scale.data["Slamf6",],
         method = "spearman")

cor(exp_agg_CD8_PE2@assays$integrated@scale.data["Tcf7",],exp_agg_CD8_PE2@assays$integrated@scale.data["Slamf6",])
cor.test(exp_agg_CD8_PE2@assays$integrated@scale.data["Tcf7",],
         exp_agg_CD8_PE2@assays$integrated@scale.data["Slamf6",],
         method = "spearman")

# Idents(exp_agg_CD8) <- "AR_bin"
# DE_AR_CD8 <- FindMarkers(exp_agg_CD8,
#                             ident.1 = "AR-",
#                             ident.2 = "AR+",
#                             logfc.threshold = 0,
#                             min.cells.group = 1,
#                             min.cells.feature = 1,
#                             min.pct = 0)
# write.csv(DE_AR_CD8, file = "figures/DE_AR/DE_AR_CD8T.csv")
DE_AR_CD8 <- read.csv("figures/DE_AR/DE_AR_CD8T.csv") %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Up-regulated","Down-regulated"),
         sig = ifelse(p_val_adj < 0.05, TRUE, FALSE)) %>%
  filter(X != "Ar")
DE_AR_CD8_sig <- DE_AR_CD8 %>%
  filter(p_val_adj < 0.05)
write.csv(DE_AR_CD8_sig, file = "figures/DE_AR/DE_AR_CD8T_sig.csv")
top_DE_AR_CD8_up <- DE_AR_CD8 %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_DE_AR_CD8_down <- DE_AR_CD8 %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
label_genes_DE_AR_CD8 <- c(top_DE_AR_CD8_up,top_DE_AR_CD8_down)
v_AR_CD8 <- ggplot(data = DE_AR_CD8,
                      aes(x = (avg_log2FC),
                          y = -log10(p_val_adj + .Machine$double.xmin),
                          color = sig)) +
  geom_point(alpha = 0.6) +
  # annotate("label", x = 1.8, y = 18, label = "Up in AR-")+
  # annotate("label", x = -1.8, y = 18, label = "Up in AR+") + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) +
  geom_label_repel(aes(label = ifelse(X %in% label_genes_DE_AR_CD8,as.character(X),""),
                       color = pol)) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(expand = c(0.01,0.01)) +
  #scale_x_continuous(limits = c(-2.5,2.5)) + 
  xlab("Average log-fold change") +
  ylab("-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0) 
ggsave(v_AR_CD8,
       filename = "figures/AR_CD8_volcano.pdf",
       height = 6,
       width = 8)



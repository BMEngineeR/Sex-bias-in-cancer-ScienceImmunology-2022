library(Seurat)

# load data & create Seurat objects
# focus on just BBN mouse data

# control - 10 weeks 
control_10 <- Read10X_h5(filename = "data/Murine_BBN/control_10 weeks_raw_feature_bc_matrix.h5")
exp_con_10 <- CreateSeuratObject(control_10,
                                 project = "CON_10",
                                 min.cells = 5)
exp_con_10 <- subset(exp_con_10, subset = nFeature_RNA > 500)
exp_con_10 <- NormalizeData(exp_con_10, verbose = FALSE)
exp_con_10 <- FindVariableFeatures(exp_con_10, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
exp_con_10$group <- "control"
exp_con_10$week <- "10"
exp_con_10$replicate <- "1"

# control - 23 weeks 
control_23 <- Read10X_h5(filename = "data/Murine_BBN/control_23 weeks_raw_feature_bc_matrix.h5")
exp_con_23 <- CreateSeuratObject(control_23,
                                 project = "CON_23",
                                 min.cells = 5)
exp_con_23 <- subset(exp_con_23, subset = nFeature_RNA > 500)
exp_con_23 <- NormalizeData(exp_con_23, verbose = FALSE)
exp_con_23 <- FindVariableFeatures(exp_con_23, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
exp_con_23$group <- "control"
exp_con_23$week <- "23"
exp_con_23$replicate <- "1"

# BBN - 10.5 weeks
BBN_10.5 <- Read10X_h5(filename = "data/Murine_BBN/BBN_10.5 weeks_raw_feature_bc_matrix.h5")
exp_BBN_10.5 <- CreateSeuratObject(BBN_10.5,
                                   project = "BBN_10.5",
                                   min.cells = 5)
exp_BBN_10.5 <- subset(exp_BBN_10.5, subset = nFeature_RNA > 500)
exp_BBN_10.5 <- NormalizeData(exp_BBN_10.5, verbose = FALSE)
exp_BBN_10.5 <- FindVariableFeatures(exp_BBN_10.5, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)
exp_BBN_10.5$group <- "BBN"
exp_BBN_10.5$week <- "10.5"
exp_BBN_10.5$replicate <- "1"

# BBN - 17.5 weeks #1
BBN_17.5_1 <- Read10X_h5(filename = "data/Murine_BBN/BBN_17.5 week_tumor #1_raw_feature_bc_matrix.h5")
exp_BBN_17.5_1 <- CreateSeuratObject(BBN_17.5_1,
                                     project = "BBN_17.5_1",
                                     min.cells = 5)
exp_BBN_17.5_1 <- subset(exp_BBN_17.5_1, subset = nFeature_RNA > 500)
exp_BBN_17.5_1 <- NormalizeData(exp_BBN_17.5_1, verbose = FALSE)
exp_BBN_17.5_1 <- FindVariableFeatures(exp_BBN_17.5_1, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
exp_BBN_17.5_1$group <- "BBN"
exp_BBN_17.5_1$week <- "17.5"
exp_BBN_17.5_1$replicate <- "1"

# BBN - 17.5 weeks #2
BBN_17.5_2 <- Read10X_h5(filename = "data/Murine_BBN/BBN_17.5 week_tumor #2_raw_feature_bc_matrix.h5")
exp_BBN_17.5_2 <- CreateSeuratObject(BBN_17.5_2,
                                     project = "BBN_17.5_2",
                                     min.cells = 5)
exp_BBN_17.5_2 <- subset(exp_BBN_17.5_1, subset = nFeature_RNA > 500)
exp_BBN_17.5_2 <- NormalizeData(exp_BBN_17.5_2, verbose = FALSE)
exp_BBN_17.5_2 <- FindVariableFeatures(exp_BBN_17.5_2, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
exp_BBN_17.5_2$group <- "BBN"
exp_BBN_17.5_2$week <- "17.5"
exp_BBN_17.5_2$replicate <- "2"

# BBN - 23 weeks #1 
BBN_23_1 <- Read10X_h5(filename = "data/Murine_BBN/BBN_23 weeks_tumor #1_raw_feature_bc_matrix.h5")
exp_BBN_23_1 <- CreateSeuratObject(BBN_23_1,
                                     project = "BBN_23_1",
                                     min.cells = 5)
exp_BBN_23_1 <- subset(exp_BBN_23_1, subset = nFeature_RNA > 500)
exp_BBN_23_1 <- NormalizeData(exp_BBN_23_1, verbose = FALSE)
exp_BBN_23_1 <- FindVariableFeatures(exp_BBN_23_1, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
exp_BBN_23_1$group <- "BBN"
exp_BBN_23_1$week <- "23"
exp_BBN_23_1$replicate <- "1"

# BBN - 23 weeks #2
BBN_23_2 <- Read10X_h5(filename = "data/Murine_BBN/BBN_23 weeks_tumor #2_raw_feature_bc_matrix.h5")
exp_BBN_23_2 <- CreateSeuratObject(BBN_23_2,
                                   project = "BBN_23_2",
                                   min.cells = 5)
exp_BBN_23_2 <- subset(exp_BBN_23_2, subset = nFeature_RNA > 500)
exp_BBN_23_2 <- NormalizeData(exp_BBN_23_2, verbose = FALSE)
exp_BBN_23_2 <- FindVariableFeatures(exp_BBN_23_2, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)
exp_BBN_23_2$group <- "BBN"
exp_BBN_23_2$week <- "23"
exp_BBN_23_2$replicate <- "2"

exp <- list(exp_con_10,
            exp_con_23,
            exp_BBN_10.5,
            exp_BBN_17.5_1,
            exp_BBN_17.5_2,
            exp_BBN_23_1,
            exp_BBN_23_2)

# Integrate the Seurat objects
exp_anchors <- FindIntegrationAnchors( object.list = exp )
exp_agg <- IntegrateData( anchorset = exp_anchors )

save(exp_agg,
     file = "data/Murine_BBN/integrated_seurat.RData")

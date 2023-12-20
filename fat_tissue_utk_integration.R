library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(patchwork)
set.seed(412)

fatRES2.data <- Read10X(data.dir = "fattissue_data/RES_raw_feature_bc_matrix")#sparse-matrix representation

fatCONTR2.data <- Read10X(data.dir = "fattissue_data/CONTROL_raw_feature_bc_matrix")#sparse-matrix representation

fatRES2 <- CreateSeuratObject(counts = fatRES2.data, project = "fat_res", min.cells = 3, min.features = 200)

fatCONTR2 <- CreateSeuratObject(counts = fatCONTR2.data, project = "fat_control", min.cells = 3, min.features = 200)
#### unique the cell names for combined data
cell_names_res <- colnames(fatRES2)
cell_names_con <- colnames(fatCONTR2)
head(cell_names_res)
head(cell_names_con)
combined_colnames <- make.unique(c(cell_names_res, cell_names_con))
colnames(fatRES2) <- combined_colnames[1:length(cell_names_res)]
colnames(fatCONTR2) <- combined_colnames[(length(cell_names_res) + 1):(length(cell_names_res) +length(cell_names_con))]

fatRES2[["percent.mt"]] <- PercentageFeatureSet(fatRES2, pattern = "^mt-")
fatCONTR2[["percent.mt"]] <- PercentageFeatureSet(fatCONTR2, pattern = "^mt-")





fatRES2 <- subset(fatRES2, subset = nFeature_RNA>200 & nFeature_RNA < 2000 & percent.mt <10)

fatCONTR2 <- subset(fatCONTR2, subset = nFeature_RNA>200 & nFeature_RNA < 2000 & percent.mt <10)


#normalize data
fatRES2 <- NormalizeData(fatRES2, normalization.method = "LogNormalize", scale.factor = 10000)
#feature selection(find highly variable features)
fatRES2 <- FindVariableFeatures(fatRES2, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10_res2 <- head(VariableFeatures(fatRES2), 10)

fatCONTR2 <- NormalizeData(fatCONTR2, normalization.method = "LogNormalize", scale.factor = 10000)
#feature selection(find highly variable features)
fatCONTR2 <- FindVariableFeatures(fatCONTR2, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10_con2 <- head(VariableFeatures(fatCONTR2), 10)

#scale to mean 0 and variance 1
all.genes.res2 <- rownames(fatRES2)
fatRES2<-ScaleData(fatRES2, features=all.genes.res2, verbose=F)#normalization
#regress out percentage number of mitochondrial genes (regression Correction)
fatRES2 <-ScaleData(fatRES2, 
                    #features=all.genes.res2,
                    vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)

#scale to mean 0 and variance 1
all.genes.con2 <- rownames(fatCONTR2)
fatCONTR2<-ScaleData(fatCONTR2, features=all.genes.con2, verbose=F)
#regress out percentage number of mitochondrial genes
fatCONTR2 <-ScaleData(fatCONTR2, 
                      #features=all.genes.con2, 
                      vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)

# Label different dataset
fatRES2$label <- 'RESVERATROL'
fatCONTR2$label <- 'CONTROL'



# Identify integration anchors
anchors <- FindIntegrationAnchors(object.list = list(fatRES2, fatCONTR2), dims = 1:20)
# Integrate data
seurat_combined <- IntegrateData(anchorset = anchors, dims = 1:20)

# scale after integration
all_genes_combine <- rownames(seurat_combined)
seurat_combined <- ScaleData(seurat_combined, features=all_genes_combine,verbose = FALSE)
seurat_combined <-ScaleData(seurat_combined, 
                            #features=all_genes_combine,
                            vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)







# PCA
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))
#visualize reduced dimensions
VizDimLoadings(seurat_combined, dims = 1:10, reduction = "pca")
#dimension heat map(first dimension- primary sources for heterogen.)
#cells and features ordered by PCA scores
DimHeatmap(seurat_combined, dims = 1:20, cells = 500, balanced = TRUE)
#Determine "dimensionality" of dataset (paper used 12)
ElbowPlot(seurat_combined, ndims=20, reduction="pca")


# Clustering
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:15)
seurat_combined <- FindClusters(seurat_combined
                                ,resolution = 0.5 # resolution will impact how many groups resulted
                                )
# UMAP or t-SNE
seurat_combined <- RunUMAP(seurat_combined, dims = 1:15)
seurat_combined <- RunTSNE(seurat_combined, dims = 1:15)

DimPlot(seurat_combined, reduction = "umap",group.by="label")
TSNEPlot(seurat_combined, reduction = "tsne",group.by="label")


p1 <-DimPlot(seurat_combined, reduction = "umap", label=T)
p2 <-TSNEPlot(seurat_combined, reduction="tsne", label=T)
p1_2 <- DimPlot(seurat_combined, reduction = "umap", label=T,split.by = "label")
p2_2 <-TSNEPlot(seurat_combined, reduction="tsne", label=T,split.by = "label")
p1+p1_2
p2+p2_2



saveRDS(seurat_combined, file = "./fattissue_data/results/filtering/seurat_combined10.rds")
seurat_combined <- readRDS("./fattissue_data/results/filtering/seurat_combined10.rds")

# Find differentially expressed genes
fatcombine.markers <- FindAllMarkers(seurat_combined, only.pos = TRUE, min.pct = 0.25,
                                 logfc.threshold = 0.25)
fatcombine.MARKERS<-as.data.frame(fatcombine.markers %>%
                                group_by(cluster) %>%
                                slice_max(n = 20, order_by = avg_log2FC)
)
write.csv(fatcombine.MARKERS, file = "./fattissue_data/results/filtering/fatcombine_MARKERS.csv", row.names = FALSE)

# 
# # Feature plot for specific genes
# FeaturePlot(seurat_combined, features = c("gene1", "gene2"))
# '''
# 
# 








###### analyze DCN and UNP genes
indices_res <- grep("dcn", all.genes.res2,ignore.case = TRUE)
indices_con <- grep("dcn", all.genes.con2,ignore.case = TRUE)
indices_res_ucp <- grep("ucp", all.genes.res2,ignore.case = TRUE)
indices_con_ucp <- grep("ucp", all.genes.con2,ignore.case = TRUE)

dcn_elements_res <- all.genes.res2[indices_res]
dcn_elements_con <- all.genes.con2[indices_con]
ucp_elements_res <- all.genes.res2[indices_res_ucp]
ucp_elements_con <- all.genes.con2[indices_con_ucp]


FeaturePlot(seurat_combined, features = dcn_elements_res, min.cutoff = "q9")
FeaturePlot(seurat_combined, features = ucp_elements_res, min.cutoff = "q9")

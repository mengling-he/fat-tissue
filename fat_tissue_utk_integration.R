library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(patchwork)
set.seed(412)

# read in the rds file--------------------------
merged_fat <- readRDS("./results/merged_fat.rds")


# QC & filtering -----------------------

# calculate mitochondrial percentage
merged_fat$mitoPercent <- PercentageFeatureSet(merged_fat, pattern='^mt-')

# Visualize QC metrics as a violin plot
Idents(merged_fat) <- "fat"
VlnPlot(merged_fat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent")
        #, ncol = 3
        ,group.by = NULL
)
Idents(merged_fat) <- merged_fat@meta.data$orig.ident
plot1 <- FeatureScatter(merged_fat, feature1 = "nCount_RNA", feature2 = "mitoPercent")
plot2 <- FeatureScatter(merged_fat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+geom_smooth(method = 'lm')
plot1 + plot2


merged_fat2 <- subset(merged_fat, subset = nFeature_RNA>200 & nFeature_RNA < 5000 & mitoPercent <5)
# merged_fat2
# merged_fat


# perform standard workflow steps to figure out if we see any batch effects --------
#normalize data
merged_fat2 <- NormalizeData(merged_fat2, normalization.method = "LogNormalize", scale.factor = 10000)
#feature selection(find highly variable features)
merged_fat2 <- FindVariableFeatures(merged_fat2, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_fat2), 10)

# scale to mean 0 and variance 1
all.genes <- rownames(merged_fat2)
merged_fat2<-ScaleData(merged_fat2, features=all.genes, verbose=F)#normalization
# pca
merged_fat2 <- RunPCA(object = merged_fat2)
ElbowPlot(merged_fat2)

# clustering
merged_fat2 <- FindNeighbors(object = merged_fat2, dims = 1:20
                             ,reduction = "pca")
merged_fat2 <- FindClusters(object = merged_fat2,resolution=0.5
                            , cluster.name = "unintegrated_clusters")
merged_fat2 <- RunUMAP(merged_fat2, dims = 1:20, reduction = "pca"
                       , reduction.name = "umap.unintegrated")
plot1 <- DimPlot(merged_fat2, reduction = "umap.unintegrated"
                 ,group.by = c("orig.ident", "unintegrated_clusters")
                 )
plot1
saveRDS(merged_fat2, file = "./results/merged_fat2.rds")



# Perform streamlined (rpca) integrative analysis--------------------
fat_integration <- IntegrateLayers(
  object = merged_fat2, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.rpca',
  verbose = FALSE)
# clustering
fat_integration <- FindNeighbors(fat_integration, reduction = "integrated.rpca", dims = 1:20)
fat_integration <- FindClusters(fat_integration
                                , resolution =c(0.1,0.2,0.3,0.5,0.7)
                                , cluster.name = c("rpca_clusters_0.1","rpca_clusters_0.2","rpca_clusters_0.3","rpca_clusters_0.5","rpca_clusters_0.7"))
fat_integration <- RunUMAP(fat_integration, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")


# rpca integration -  umap by orig.ident
p1 <- DimPlot(
  fat_integration,
  reduction = "umap.rpca",
  group.by = c("orig.ident"),
  combine = FALSE, label.size = 2,label = T
)

# clustering result on umap vs. rpca umap 
DimPlot(fat_integration, reduction = "umap.rpca"
        , group.by = c("orig.ident","rpca_clusters_0.1","rpca_clusters_0.2"
                       ,"rpca_clusters_0.3","rpca_clusters_0.5","rpca_clusters_0.7")
        ,label = T)



DimPlot(fat_integration, reduction = "umap.unintegrated"
        , group.by = c("orig.ident","rpca_clusters_0.1","rpca_clusters_0.2"
                       ,"rpca_clusters_0.3","rpca_clusters_0.5","rpca_clusters_0.7")
        ,label = T)

DimPlot(fat_integration, reduction = "umap.rpca"
        , group.by = c("rpca_clusters_0.2")
        ,label = T)|DimPlot(fat_integration, reduction = "umap.unintegrated"
                            , group.by = c("rpca_clusters_0.2")
                            ,label = T)


# Perform streamlined (cca) integrative analysis--------------------
fat_integration2 <- IntegrateLayers(
  object = merged_fat2, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.cca',
  verbose = FALSE)
# clustering
fat_integration2 <- FindNeighbors(fat_integration2, reduction = "integrated.cca", dims = 1:20)
fat_integration2 <- FindClusters(fat_integration2
                                , resolution =c(0.1,0.2,0.3,0.5,0.7)
                                , cluster.name = c("cca_clusters_0.1","cca_clusters_0.2","cca_clusters_0.3","cca_clusters_0.5","cca_clusters_0.7"))
fat_integration2 <- RunUMAP(fat_integration2, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")


DimPlot(fat_integration2, reduction = "umap.cca"
        , group.by = c("cca_clusters_0.1","cca_clusters_0.2"
                       ,"cca_clusters_0.3","cca_clusters_0.5")
        ,label = T)

DimPlot(fat_integration, reduction = "umap.rpca"
        , group.by = c("rpca_clusters_0.2")
        ,label = T)|DimPlot(fat_integration2, reduction = "umap.cca"
                            , group.by = c("cca_clusters_0.3")
                            ,label = T)
DimPlot(fat_integration2, reduction = "umap.cca"
        , group.by = c("cca_clusters_0.3")
        ,label = T,split.by = "orig.ident")






# analyze DCN and UNP genes--------------------
indices_dcn <- grep("dcn", all.genes,ignore.case = TRUE)
indices_ucp <- grep("ucp", all.genes,ignore.case = TRUE)

dcn_elements <- all.genes[indices_dcn]
ucp_elements <- all.genes[indices_ucp]

# DCN in cluster 1,3,4, espeacially in cluster 1 res group
VlnPlot(fat_integration2, features=dcn_elements,split.by = "orig.ident",combine=FALSE)
DimPlot(fat_integration2, reduction = "umap.cca"
        , group.by = c("cca_clusters_0.3")
        ,label = T)|FeaturePlot(fat_integration2,reduction =  "umap.cca", features = dcn_elements, min.cutoff = "q9")
FeaturePlot(fat_integration2,reduction =  "umap.cca",features = dcn_elements,split.by = "orig.ident" ,min.cutoff = "q9")

# Ucp3 in cluster 5(control)
# Ucp2 in cluster 0(res),2(res),9(both)
# Ucp1 low in all clusters
VlnPlot(fat_integration2, features=ucp_elements,combine = FALSE)
VlnPlot(fat_integration2, features=ucp_elements,combine = FALSE,split.by = "orig.ident")
FeaturePlot(fat_integration2,reduction =  "umap.cca",features = "Ucp1",split.by = "orig.ident" ,min.cutoff = "q9")
# # Warning message:
# All cells have the same value (1.3652867355258) of “Ucp1” 






# findConserved markers (integration= cca ; resolution=0.3)--------------------
fat_integration2 <- JoinLayers(fat_integration2)
#saveRDS(fat_integration2, file = "./results/fat_integration2.rds")

DimPlot(fat_integration2, reduction = "umap.cca"
        , group.by = c("orig.ident","cca_clusters_0.3")
        ,label = T)

Idents(fat_integration2) <- fat_integration2@meta.data$cca_clusters_0.3

n_cluster <- length(levels(Idents(fat_integration2)))

marker = list()
for (i in 1:n_cluster){
  marker[[i]] <- FindConservedMarkers(fat_integration2, ident.1 = i-1, grouping.var = "orig.ident", verbose = FALSE)
}


head(marker[[1]],5)

N=20
#cluster <- rep(0:(n_cluster-1), each = N)
marker_list_20 <- lapply(marker, function(df) head(df,N))
name_of_markers <- lapply(marker_list_20, function(df) rownames(df))
name_of_markers <- unlist(name_of_markers)
marker.df <- bind_rows(marker_list_20, .id = 'cluster')
marker.df$cluster <- as.integer(marker.df$cluster)-1
marker.df$cluster <- as.character(marker.df$cluster)
marker.df$marker <- name_of_markers
head(marker.df)
write.csv(marker.df, "./results/genemarker_updated.csv", row.names=FALSE)

# check duplicated genemarker
marker.df <- read.csv("./results/genemarker_updated.csv")
d <- as.data.frame(table(marker.df$marker))

duplicated_rows <- marker.df[duplicated(marker.df$marker) | duplicated(marker.df$marker, fromLast = TRUE), ]
d2 <- table(duplicated_rows$marker,duplicated_rows$cluster)
# custer 1 and cluster 9 has many 10 samee gene markers, and clustrer 3 and 4 has 1 same gene marker(Slit3)
duplicated_rows$marker
VlnPlot(fat_integration2, features='Opcml',combine = FALSE)


# clutsering label----------------
cluster_bycell <- fat_integration2@meta.data[c("orig.ident","cca_clusters_0.3")]
tab = table(cluster_bycell)
tab <- cbind(tab, Total = rowSums(tab))

write.csv(cluster_bycell, "./results/cluster_bycell.csv", row.names=TRUE)





# findMarkers between conditions ---------------------
fat_integration2$cluster.sample <- paste0(fat_integration2$cca_clusters_0.3,'.', fat_integration2$orig.ident)
View(fat_integration2@meta.data)
Idents(fat_integration2) <- fat_integration2$cluster.sample
DimPlot(fat_integration2, reduction = 'umap.cca', label = TRUE)


interferon_list = list()
for (i in 1:n_cluster){
  ident1 <- paste0(i-1,'.','fat_res')
  ident2 <- paste0(i-1,'.','fat_control')
  interferon_list[[i]] <- FindMarkers(fat_integration2, ident.1 = ident1, ident.2 = ident2)
}

interferon_list_20 <- lapply(interferon_list, function(df) head(df,20))
interferon_list_20[[1]]
rowname_of_interferon <- lapply(interferon_list_20, function(df) rownames(df))
rowname_of_interferon <- unlist(rowname_of_interferon)
interferon_20 <- bind_rows(interferon_list_20, .id = "cluster")
interferon_20$cluster <- as.integer(interferon_20$cluster)-1
interferon_20$cluster <- as.character(interferon_20$cluster)
interferon_20$marker <- rowname_of_interferon
head(interferon_20)

write.csv(interferon_20, "./results/genemarker_bysample.csv", row.names=FALSE)
# compare gene markers (cluster 0)
FeaturePlot(fat_integration2, features = c('F13a1', 'Mctp1', 'Mrc1'), split.by = 'orig.ident', min.cutoff = 'q10',reduction = 'umap.cca')
FeaturePlot(fat_integration2, features = c('Fabp4', 'Cfd', 'Fkbp5'), split.by = 'orig.ident', min.cutoff = 'q10',reduction = 'umap.cca')

























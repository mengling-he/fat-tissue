library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(patchwork)
set.seed(412)


# read data
fatRES2.data <- Read10X(data.dir = "./RES_raw_feature_bc_matrix")#sparse-matrix representation
fatCONTR2.data <- Read10X(data.dir = "./CONTROL_raw_feature_bc_matrix")#sparse-matrix representation

# create seurat
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

PercentageFeatureSet(fatRES2, pattern='^mt-')



# merge Seurat -----------------------
merged_fat <- merge(x = fatRES2, y = fatCONTR2)
merged_fat
View(merged_fat@meta.data)

saveRDS(merged_fat, file = "./results/merged_fat.rds")


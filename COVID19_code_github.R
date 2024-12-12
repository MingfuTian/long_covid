# Figure 1


#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## Work name
work_name <- "COVID19"
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## input data
covid_21_035_count <- Read10X_h5(filename = "InputData/21_035/21_035_filtered_feature_bc_matrix.h5")
covid_21_035 <- CreateSeuratObject(counts = covid_21_035_count,assay = 'RNA',min.cells= 10, min.features = 500)
covid_21_035$orig.ident  <- "covid1"
covid_21_035$sample_name <- "covid_21_035"

covid_BB_22001_count <- Read10X_h5(filename = "InputData/BB_22001/BB_22001_filtered_feature_bc_matrix.h5")
covid_BB_22001 <- CreateSeuratObject(counts = covid_BB_22001_count,assay = 'RNA',min.cells= 10, min.features = 500)
covid_BB_22001$orig.ident <- "covid2"
covid_BB_22001$sample_name <- "covid_BB_22001"

covid_LIBB_0093_count <- Read10X_h5(filename = "InputData/LIBB_0093/LIBB_0093_filtered_feature_bc_matrix.h5")
covid_LIBB_0093 <- CreateSeuratObject(counts = covid_LIBB_0093_count,assay = 'RNA',min.cells= 10, min.features = 500)
covid_LIBB_0093$orig.ident <- "covid3"
covid_LIBB_0093$sample_name <- "covid_LIBB_0093"

covid_CC003_count <- Read10X_h5(filename = "InputData/CC003/CC003_filtered_feature_bc_matrix.h5")
covid_CC003 <- CreateSeuratObject(counts = covid_CC003_count,assay = 'RNA',min.cells= 10, min.features = 500)
covid_CC003$orig.ident <- "ctrl1"
covid_CC003$sample_name <- "covid_CC003"

covid_CC005_19_count <- Read10X_h5(filename = "InputData/CC005_19/CC005_19_filtered_feature_bc_matrix.h5")
covid_CC005_19 <- CreateSeuratObject(counts = covid_CC005_19_count,assay = 'RNA',min.cells= 10, min.features = 500)
covid_CC005_19$orig.ident <- "ctrl2"
covid_CC005_19$sample_name <- "covid_CC005_19"

## merge sample
covid19 <- merge(covid_21_035, y = c(covid_BB_22001, covid_LIBB_0093, covid_CC003,covid_CC005_19), add.cell.ids = c("covid1", "covid2", "covid3", "ctrl1","ctrl2"))
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## QC
covid19[["Percent_MT_Gene"]] <- PercentageFeatureSet(covid19, pattern = "^MT-")
VlnPlot(covid19, features = c("nFeature_RNA", "nCount_RNA", "Percent_MT_Gene"), ncol = 3)
covid19 <- subset(covid19, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 200 & nCount_RNA < 25000 & Percent_MT_Gene < 10)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## cell cluster
covid19 <- JoinLayers(covid19)
covid19[["RNA"]] <- split(covid19[["RNA"]], f = covid19$orig.ident)
covid19 <- NormalizeData(covid19, normalization.method = "LogNormalize", scale.factor = 10000)
covid19 <- FindVariableFeatures(covid19, selection.method = "vst", nfeatures = 2000)
covid19 <- ScaleData(covid19, split.by = "orig.ident"	, features = rownames(covid19))
covid19 <- RunPCA(covid19, features = VariableFeatures(object = covid19))
covid19 <- FindNeighbors(covid19, dims = 1:30, reduction = "pca")
covid19 <- FindClusters(covid19, resolution = 2)
covid19 <- RunUMAP(covid19, dims = 1:30, reduction = "pca")

DimPlot(covid19, group.by = "seurat_clusters")
DimPlot(covid19, group.by = "orig.ident", label = T)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## harmony
library(harmony)
covid19 <- IntegrateLayers(object = covid19,method = HarmonyIntegration,orig.reduction = "pca",new.reduction = "harmony",verbose = FALSE)
covid19 <- FindNeighbors(covid19, reduction = "harmony", dims = 1:30)
covid19 <- FindClusters(covid19, resolution = 1, cluster.name = "harmony_clusters")
covid19 <- RunUMAP(covid19, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## update cluster id
covid19@meta.data <- covid19@meta.data %>% mutate(cluster_id = as.numeric(as.character(harmony_clusters)) + 1)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## cell marker
library(presto)
covid19 <- JoinLayers(covid19)
covid19_cell_marker <- FindAllMarkers(covid19, assay = "RNA", only.pos = T, logfc.threshold = 0.25,  test.use = "wilcox", min.pct = 0.1)
write.table(covid19_cell_marker, file = "covid19_cell_marker.tsv", sep = "\t", row.names = F)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## cluster anno
big_cluster_name <- c(
  `1`="B cell",
  `2`="Fibroblast",
  `3`="Epithelial cell",
  `4`="CD4+ T cell",
  `5`="Macrophage",
  `6`="CD8+ T cell",
  `7`="Fibroblast",
  `8`="NK cell",
  `9`="Epithelial cell",
  `10`="Fibroblast",
  `11`="Endothelial cell",
  `12`="Epithelial cell",
  `13`="Fibroblast",
  `14`="Epithelial cell",
  `15`="Fibroblast",
  `16`="Macrophage",
  `17`="Macrophage",
  `18`="Epithelial cell",
  `19`="Club cell",
  `20`="Epithelial cell",
  `21`="Fibroblast",
  `22`="Fibroblast"
)

## cluster color
big_cluster_pal <- c(
  "B cell"       = "#ed1299",
  "CD4+ T cell"  = "#dd27ce",
  "CD8+ T cell"  = "#273c75",
  "NK cell"      = "#c93f00",
  "Fibroblast" = "#8c7ae6",
  "Macrophage"    = "#ff523f",
  "Epithelial cell"   = "#2927c4",
  "Endothelial cell"      = "#44bd32",
  "Club cell" = "red",
  "Regulatory T cell"  = "#a1ce4c"
)

## add cluster name
covid19[['big_cluster_name']] = unname(big_cluster_name[covid19@meta.data$cluster_id])
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## cluster umap
library(scCustomize)
covid19_cluster <- UMAPPlot(covid19, reduction = "umap.harmony", group.by = "big_cluster_name", label = TRUE, label.box = T, label.color = "#000000", pt.size = 0.2, cols = big_cluster_pal)
covid19_cluster

covid19_sample <- UMAPPlot(covid19, reduction = "umap.harmony", group.by = "orig.ident", label = TRUE, label.box = T, label.color = "#000000", pt.size = 0.2)
covid19_sample
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## cell marker density
covid19_cd68 <- Plot_Density_Custom(seurat_object = covid19, features = "CD68")
covid19_cd68

covid19_MRC1 <- Plot_Density_Custom(seurat_object = covid19, features = "MRC1")
covid19_MRC1

covid19_CD163 <- Plot_Density_Custom(seurat_object = covid19, features = "CD163")
covid19_CD163

covid19_CSF1R <- Plot_Density_Custom(seurat_object = covid19, features = "CSF1R")
covid19_CSF1R

covid19_EPCAM <- Plot_Density_Custom(seurat_object = covid19, features = "EPCAM")
covid19_EPCAM

covid19_CDH1 <- Plot_Density_Custom(seurat_object = covid19, features = "CDH1")
covid19_CDH1
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
## cluster pathway
library(ReactomeGSA.data)
library(ReactomeGSA)

covid19 <- JoinLayers(covid19)
covid19_gsva <- analyse_sc_clusters(covid19, verbose = TRUE)
covid19_pathway_exp <- pathways(covid19_gsva)

relevant_pathways <- c("R-HSA-2559584", "R-HSA-2559585", "R-HSA-2559586", "R-HSA-9825895", "R-HSA-1606322","R-HSA-3134973","R-HSA-3270619","R-HSA-9732724","R-HSA-3134963","R-HSA-9758919","R-HSA-5676594")

tb_ins_covid19_pathway_exp <- covid19_pathway_exp %>% rownames_to_column(var = "id") %>% dplyr::filter(id %in% relevant_pathways)
tb_ins_covid19_pathway_exp_tmp <- tb_ins_covid19_pathway_exp[,-1] %>% column_to_rownames(var = "Name")

library(pheatmap)
pheatmap::pheatmap(tb_ins_covid19_pathway_exp_tmp, scale = "row", color = colorRampPalette(colors = c("blue","white","red"))(100))
#----------------------------------------------------------------------------------------------------

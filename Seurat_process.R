library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(future)
library(harmony)

count <- Read10X(data.dir = "F:/Single cell RNA seq/AGG_outs/AGGP14_WD_S1_2Batch/outs/count/filtered_feature_bc_matrix")
cells <- colnames(count)

#4_batch
Ctrl_S1_data1 <- count[,grepl('1$', cells)]  
Ctrl_S1_data2 <- count[,grepl('2$', cells)] 
WD_S1_data1 <- count[,grepl('3$', cells)]
WD_S1_data2 <- count[,grepl('4$', cells)]

Ctrl_S1_1_raw_cell_number <- paste('Ctrl_S1', dim(Ctrl_S1_data1)[2],' cells')    
Ctrl_S1_2_raw_cell_number <- paste('Ctrl_S1', dim(Ctrl_S1_data2)[2],' cells')
WD_S1_1_raw_cell_number <- paste('WD_S1', dim(WD_S1_data1)[2],' cells')       
WD_S1_2_raw_cell_number <- paste('WD_S1', dim(WD_S1_data2)[2],' cells')       



print(Ctrl_S1_raw_cell_number)  
print(WD_S1_raw_cell_number)


#create Seurat Object respectively
#Ctrl data
Ctrl_S1_1 <- CreateSeuratObject(counts = Ctrl_S1_data1, project = "Ctrl_S1_1",min.cells = 3, min.features = 200)
Ctrl_S1_1@meta.data$stim <- "Ctrl"
Ctrl_S1_1@meta.data$group_ID <- "Ctrl_1"

Ctrl_S1_2 <- CreateSeuratObject(counts = Ctrl_S1_data2, project = "Ctrl_S1_2",min.cells = 3, min.features = 200)
Ctrl_S1_2@meta.data$stim <- "Ctrl"
Ctrl_S1_1@meta.data$group_ID <- "Ctrl_2"

#WD data
WD_S1_1 <- CreateSeuratObject(counts = WD_S1_data1,project = "WD_S1_1", min.cells = 3, min.features = 200)
WD_S1_1@meta.data$stim <- "WD"
Ctrl_S1_1@meta.data$group_ID <- "WD_1"

WD_S1_2 <- CreateSeuratObject(counts = WD_S1_data2,project = "WD_S1_2", min.cells = 3, min.features = 200)
WD_S1_2@meta.data$stim <- "WD"
Ctrl_S1_1@meta.data$group_ID <- "WD_2"

WD.combined.raw <- merge(Ctrl_S1_1, y = c(Ctrl_S1_2,WD_S1_1,WD_S1_2),  project = "S1")


WD.combined.raw[["percent.mt"]] <- PercentageFeatureSet(WD.combined.raw, pattern = "^mt-")
WD.combined.raw[["percent.ribo"]] <- PercentageFeatureSet(WD.combined.raw, pattern = "^Rp[sl]")
WD.combined.raw[["percent.HB"]] <- PercentageFeatureSet(WD.combined.raw, pattern = "^Hb[ba]-")


VlnPlot(WD.combined.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.HB"), ncol = 5)
WD.combined.filted <- subset(WD.combined.raw, nFeature_RNA > 800 & nFeature_RNA < 9000 & nCount_RNA < 60000 & percent.mt < 5 & percent.ribo <30 & percent.HB < 0.05)
VlnPlot(WD.combined.filted, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.HB"), ncol = 5)

Cell_number <- table(WD.combined@meta.data$stim)
View(Cell_number)

#multisession
plan(sequential)
plan("multisession", workers = 18)
options(future.globals.maxSize = 1000 * 1024^2*20)

WD.combined.filted <- NormalizeData(WD.combined.filted, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","percent.HB","percent.ribo") ) %>% 
  RunPCA(verbose=TRUE)

#Run Harmony
WD.combined <- RunHarmony(WD.combined.filted,"stim",plot_convergence = T)
#Elbowplot decide dims
harmony_embeddings <- Embeddings(WD.combined, 'harmony')
harmony_embeddings[1:5, 1:5]
WD.combined <- FindNeighbors(WD.combined, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.3)
WD.combined <- RunUMAP(WD.combined,  umap.method = "umap-learn", metric = "correlation",reduction = "harmony", dims = 1:20)


p1 <- DimPlot(WD.combined, reduction = "umap", group.by = "orig.ident")+NoAxes()
p2 <- DimPlot(WD.combined, reduction = "umap", label = TRUE)+NoAxes()
plot_grid(p1, p2)

#FindConservedMarkers for cluster 0-17
sheet0 <- FindConservedMarkers(WD.combined, ident.1 = 0-17, grouping.var = "stim", verbose = TRUE)

WD.combined <- RenameIdents(WD.combined, `0` = "Microglia", `1` = "L4", `2` = "OPC", 
                            `3` = "EC", `4` = "COP", `5` = "L2/3", 
                            `6` = "Pericyte", `7` = "NFOL", `8` = "L5 IT",
                            `9` = "L6 CT", `10` = "Interneuron", `11` = "L6 IT",`12` = "OPCcyc",
                            `13` = "Astrocyte", `14` = "Pvm", `15` = "L5 NP",`16` = "Fibroblast",
                            `17` = "SMC")

WD.combined$celltype <- factor(WD.combined$celltype, levels = c("L2/3","L4",
                                                                "L5 IT","L5 NP", 
                                                                "L6 CT", "L6 IT",
                                                                "Interneuron",
                                                                "OPCcyc","OPC","COP","MFOL",
                                                                "EC","SMC","Pericyte","Fibroblast","Astrocyte","Microglia","Pvm"))

DimPlot(WD.combined, reduction = "umap", group.by="celltype",label = TRUE,cols = c("#96c3dc","#4e79a6","#abd0a7","#5fa664","#f0e2a3","#e1c548",'#ca6a6b',"#e5b5b5","#CCC9E9","#a199be","#7559A2FF",
                                                                                   "#d7ee96","#aedd2f","#91D0BE","#cdaa9f","#967568","#bac4d0","#9FA3A8"))+NoAxes()


saveRDS(WD.combined, file = "WD.combined.celltype.rds")

#Separate conditions
WD.combined$celltype.stim <- paste(Idents(WD.combined), WD.combined$stim, sep = "_")
WD.combined$celltype <- Idents(WD.combined)
Idents(WD.combined) <- "celltype.stim"

saveRDS(WD.combined, file = "WD.combined.stim.rds")

#Subset Glutamatergic neurons
WD.combined <- readRDS("WD.combined.celltype.rds")
WD.combined <- subset(WD.combined, idents = c("L2/3", "L4", "L5 IT", "L5 NP", "L6 CT", "L6 IT"))

options(future.globals.maxSize = 1000 * 1024^2*20)

WD.combined.filted <- NormalizeData(WD.combined.filted, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","percent.HB","percent.ribo") ) %>% 
  RunPCA(verbose=TRUE)

#Run Harmony
WD.combined <- RunHarmony(WD.combined.filted,"stim",plot_convergence = T)
#Elbowplot decide dims
harmony_embeddings <- Embeddings(WD.combined, 'harmony')
harmony_embeddings[1:5, 1:5]
WD.combined <- FindNeighbors(WD.combined, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.3)
WD.combined <- RunUMAP(WD.combined,  umap.method = "umap-learn", metric = "correlation",reduction = "harmony", dims = 1:20)

#FindConservedMarkers for cluster 0-6
sheet0 <- FindConservedMarkers(WD.combined, ident.1 = 0-6, grouping.var = "stim", verbose = TRUE)

WD.combined <- RenameIdents(WD.combined, `0` = "L4", `1` = "L2/3", `2` = "L6 CT", 
                            `3` = "L5 ITa", `4` = "L6 IT", `5` = "L5 ITb", `6` = "L5 NP")

WD.combined$celltype <- Idents(WD.combined)
WD.combined$celltype <- factor(WD.combined$celltype, levels = c("L2/3", "L4", "L5 ITa", "L5 ITb", "L5 NP", "L6 IT", "L6 CT"))
  
DimPlot(WD.combined, reduction = "umap", group.by="celltype",label = TRUE,cols = c("#96c3dc","#4e79a6","#abd0a7","#5fa664","#f0e2a3","#e1c548",'#ca6a6b'))

saveRDS(WD.combined, file = "WD.combined.glutamatergic.neurons.rds")

WD.combined$celltype.stim <- paste(Idents(WD.combined), WD.combined$stim, sep = "_")
WD.combined$celltype <- Idents(WD.combined)
Idents(WD.combined) <- "celltype.stim"

saveRDS(WD.combined, file = "WD.combined.glutamatergic.stim.rds")

#Augur analysis
library(Augur)

WD.combined <- readRDS("WD.combined.stim.rds")

meta_data <- WD.combined@meta.data

new_data_frame <- meta_data[, c("stim","celltype")]

colnames(new_data_frame)<-c("label","cell_type")

WD.combined@meta.data <- new_data_frame

augur = calculate_auc(WD.combined, cell_type_col = "cell_type", label_col = "label")

augur$AUC

plot_lollipop(augur)

plot_umap(
  augur,sc=WD.combined,mode = c("rank"),
  reduction = "umap",
  palette = "RdGy",
  top_n = 10,
  cell_type_col = "cell_type")

#Transfer .rds to .h5ad for SCANPY
library(sceasy)
library(reticulate)

reticulate:::find_conda()
py_available()
Sys.which("python")
loompy <- reticulate::import('loompy')

WD.combined.celltype <- readRDS("WD.combined.celltype.rds")
WD.combined.celltype[["RNA"]] <- as(WD.combined.celltype[["RNA"]], "Assay")
sceasy::convertFormat(WD.combined.celltype, from="seurat", to="anndata", outFile='/WD.combined.celltype.h5ad')

WD.combined.Glu.stim <- readRDS("WD.combined.glutamatergic.stim.rds")
WD.combined.Glu.stim[["RNA"]] <- as(WD.combined.Glu.stim[["RNA"]], "Assay")
sceasy::convertFormat(WD.combined.Glu.stim, from="seurat", to="anndata", outFile='/WD.combined.Glu.stim.h5ad')

WD.combined <- readRDS("WD.combined.stim.rds")
WD.combined.L234 <- subset(WD.combined, idents = c("L4_WD","L4_Ctrl",
                                              "L2/3_WD","L2/3_Ctrl"))
WD.combined.L234[["RNA"]] <- as(WD.combined.L234[["RNA"]], "Assay")
sceasy::convertFormat(WD.combined.L234, from="seurat", to="anndata", outFile='WD_neuron_L234.h5ad')

WD.combined.L56 <- subset(WD.combined, idents = c("L5 ITa_WD","L5 ITa_Ctrl",
                                                  "L5 ITb_WD","L5 ITb_Ctrl",
                                                  "L5 NP_WD","L5 NP_Ctrl",
                                                  "L6 CT_WD","L6 CT_Ctrl",
                                                  "L6 IT_WD","L6 IT_Ctrl"))
WD.combined.L56[["RNA"]] <- as(WD.combined.L56[["RNA"]], "Assay")
sceasy::convertFormat(WD.combined.L56, from="seurat", to="anndata", outFile='WD_neuron_L56.h5ad')

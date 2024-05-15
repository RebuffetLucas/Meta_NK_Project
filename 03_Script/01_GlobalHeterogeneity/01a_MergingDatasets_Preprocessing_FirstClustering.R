#Merging All the datasets together and running the different batch corrections and Obtaining the first clustering

#OpenAll

#ForV2 samples only:
PBMC = readRDS(PATH_MERGED_DOWNSAMPLE_Obj)
PBMC= subset(PBMC, subset= Chemistry == "V2" )


#Merge them
Merged_Seurat= NoIPH
Merged_Seurat$orig.ident = droplevels(Merged_Seurat$orig.ident)


#Free Memory
gc()


#Rescale with MultiBatchNorm (take into account the sequencing depth of each Dataset
# !!! This step replace NormalizeDataStep !!!
sce_Merged= as.SingleCellExperiment(Merged_Seurat)
sce_Merged_Rescaled = multiBatchNorm(sce_Merged,  batch = sce_Merged$orig.ident)
Merged_Seurat_Rescaled= as.Seurat(sce_Merged_Rescaled)


#Find VariableFeatures across samples
#Variable on orig.ident
seurat_list_Rescaled = SplitObject(Merged_Seurat_Rescaled, split.by = "orig.ident")
seurat_list_Rescaled <- lapply(seurat_list_Rescaled, function(x) FindVariableFeatures(x, nfeatures= 5000, selection.method = "vst"))
VARIABLE_FEATURES = SelectIntegrationFeatures(  seurat_list_Rescaled,  nfeatures = 2000,  assay = NULL,  verbose = TRUE)


#Free Memory
rm(Merged_Seurat , sce_Merged , sce_Merged_Rescaled , seurat_list_Rescaled )
gc()


#Run clustering

All_genes= rownames(Merged_Seurat_Rescaled)
Merged_Seurat_Rescaled=ScaleData(Merged_Seurat_Rescaled, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )

Merged_Seurat_Rescaled=RunPCA(Merged_Seurat_Rescaled, features= VARIABLE_FEATURES)


#Run batch Correction
Merged_Seurat_Rescaled=Merged_Seurat_Rescaled %>%                  
  RunHarmony(c("orig.ident"), plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:30 ) 

#Run UMAP
Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30, reduction.name = "UMAP", reduction.key = "UMAP")


#Run clustering
Merged_Seurat_Rescaled <- FindNeighbors(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30)
Merged_Seurat_Rescaled <- FindClusters(Merged_Seurat_Rescaled, resolution = 0.7, random.seed = SEED )

DimPlot(Merged_Seurat_Rescaled, label= TRUE) #Have a first look


#RemoveDying and Prolif
Merged_Seurat_Rescaled = subset(Merged_Seurat_Rescaled , idents= c("0", "1", "2", "3", "4", "5", "6", "8"))




#Rerun the pipeline on the curated datasets
Merged_Seurat_Rescaled=RunPCA(Merged_Seurat_Rescaled, features= VARIABLE_FEATURES)

DimPlot(Merged_Seurat_Rescaled, reduction = "pca", group.by = "Dataset")
VizDimLoadings(Merged_Seurat_Rescaled, reduction = "pca")

#Free Memory
gc()

Merged_Seurat_Rescaled=Merged_Seurat_Rescaled %>%                  
  RunHarmony(c("orig.ident"), plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:30 ) 

Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30, reduction.name = "UMAP", reduction.key = "UMAP")







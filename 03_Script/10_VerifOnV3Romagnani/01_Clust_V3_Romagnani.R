#Merging All the datasets together and running the different batch corrections and Obtaining the first clustering

#OpenAll

#ForV2:
NoIPH = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/00_RawData/Seurat_Obj_other_samples/Sent_By_Janine/merged_and_downsampled.rds")
NoIPH= subset(NoIPH, subset= Chemistry == "V3" )


#Merge them
Merged_Seurat= NoIPH
Merged_Seurat$orig.ident = droplevels(Merged_Seurat$orig.ident)


#Free Memory
gc()


#Rescale with MultiBatchNorm (take into account the sequencing depth of each Dataset)
# !!! This step replace NormalizeDataStep !!!
sce_Merged= as.SingleCellExperiment(Merged_Seurat)
sce_Merged_Rescaled = multiBatchNorm(sce_Merged,  batch = sce_Merged$orig.ident)
Merged_Seurat_Rescaled= as.Seurat(sce_Merged_Rescaled)

#Only LogNorm
#Merged_Seurat_Rescaled = NormalizeData(Merged_Seurat) #Do not run this step if MultiBatchNorm was used


#Find VariableFeatures across samples
#Variable on orig.ident
seurat_list_Rescaled = SplitObject(Merged_Seurat_Rescaled, split.by = "orig.ident")
seurat_list_Rescaled <- lapply(seurat_list_Rescaled, function(x) FindVariableFeatures(x, nfeatures= 5000, selection.method = "vst"))
VARIABLE_FEATURES = SelectIntegrationFeatures(  seurat_list_Rescaled,  nfeatures = 2000,  assay = NULL,  verbose = TRUE)

#Select only the variable features shared by all samples
VARIABLE_FEATURES2 = lapply(seurat_list_Rescaled, function(x) VariableFeatures(x))
VARIABLE_FEATURES2 = Reduce(intersect, VARIABLE_FEATURES2)

#Variable = Those that appear in more than 18 Datasets with 5000
VARIABLE_FEATURES5 = lapply(seurat_list_Rescaled, function(x) VariableFeatures(x))
VARIABLE_FEATURES5 = unlist(VARIABLE_FEATURES5)
names(VARIABLE_FEATURES5)=NULL
tableau = table(VARIABLE_FEATURES5)
tableau = sort(tableau, decreasing = TRUE)
table(tableau)
VARIABLE_FEATURES5 = names(tableau[tableau>=19])


#Variable = Those that appear in more than 16 Datasets with 2000
seurat_list_Rescaled = SplitObject(Merged_Seurat_Rescaled, split.by = "orig.ident")
seurat_list_Rescaled <- lapply(seurat_list_Rescaled, function(x) FindVariableFeatures(x, nfeatures= 2000, selection.method = "vst"))
VARIABLE_FEATURES6 = lapply(seurat_list_Rescaled, function(x) VariableFeatures(x))
VARIABLE_FEATURES6 = unlist(VARIABLE_FEATURES6)
names(VARIABLE_FEATURES6)=NULL
tableau = table(VARIABLE_FEATURES6)
tableau = sort(tableau, decreasing = TRUE)
table(tableau)
VARIABLE_FEATURES6 = names(tableau[tableau>=16])





#Variable = IntegrationFeatures on Dataset
seurat_list_Rescaled = SplitObject(Merged_Seurat_Rescaled, split.by = "Dataset")
seurat_list_Rescaled <- lapply(seurat_list_Rescaled, function(x) FindVariableFeatures(x, nfeatures= 5000, selection.method = "vst"))
VARIABLE_FEATURES3 = SelectIntegrationFeatures(  seurat_list_Rescaled,  nfeatures = 2000,  assay = NULL,  verbose = TRUE)

#Variable = Intersects on Datasets
VARIABLE_FEATURES4 = lapply(seurat_list_Rescaled, function(x) VariableFeatures(x))
VARIABLE_FEATURES4 = Reduce(intersect, VARIABLE_FEATURES4)









#Free Memory
rm(Merged_Seurat , sce_Merged , sce_Merged_Rescaled , seurat_list_Rescaled )
gc()


#Find only the Variable features of the merged version
#Merged_Seurat_Rescaled = FindVariableFeatures(Merged_Seurat_Rescaled, nfeatures= 5000, selection.method = "vst")

#Vizualising the best variable features
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(Merged_Seurat_Rescaled), 100)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(Merged_Seurat_Rescaled)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2


#Run clustering

All_genes= rownames(Merged_Seurat_Rescaled)
Merged_Seurat_Rescaled=ScaleData(Merged_Seurat_Rescaled, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )

VARIABLE_FEATURES_METHOD = VARIABLE_FEATURES


Merged_Seurat_Rescaled=RunPCA(Merged_Seurat_Rescaled, features= VARIABLE_FEATURES_METHOD)

DimPlot(Merged_Seurat_Rescaled, reduction = "pca", group.by = "orig.ident")
VizDimLoadings(Merged_Seurat_Rescaled, reduction = "pca")

#Free Memory
gc()

Merged_Seurat_Rescaled=Merged_Seurat_Rescaled %>%                  
  RunHarmony(c("orig.ident"), plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:30 ) 

Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30, reduction.name = "UMAP", reduction.key = "UMAP")
#Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled,  reduction.name = "UMAP", reduction.key = "UMAP")

Merged_Seurat_Rescaled <- FindNeighbors(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30)
Merged_Seurat_Rescaled <- FindClusters(Merged_Seurat_Rescaled, resolution = 0.8, random.seed = SEED ) #0.8

DimPlot(Merged_Seurat_Rescaled, label= TRUE)

#Merged_Seurat_Rescaled$Predicted_Clust = PBMC_V3$predicted.id
#DimPlot(Merged_Seurat_Rescaled, group.by = "Predicted_Clust" ,label= TRUE)

#RemoveDying and Prolif

Merged_Seurat_Rescaled = subset(Merged_Seurat_Rescaled , idents= c("0", "1", "2", "3", "4", "5", "6", "7")) # stop at 7 included


DimPlot(Merged_Seurat_Rescaled, label= TRUE)

#Rerun clustering
Merged_Seurat_Rescaled=RunPCA(Merged_Seurat_Rescaled, features= VARIABLE_FEATURES_METHOD)

DimPlot(Merged_Seurat_Rescaled, reduction = "pca", group.by = "Dataset")
VizDimLoadings(Merged_Seurat_Rescaled, reduction = "pca")

#Free Memory
gc()

Merged_Seurat_Rescaled=Merged_Seurat_Rescaled %>%                  
  RunHarmony(c("orig.ident"), plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:30 ) 

Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30, reduction.name = "UMAP", reduction.key = "UMAP")
#Merged_Seurat_Rescaled <- RunUMAP(Merged_Seurat_Rescaled,  reduction.name = "UMAP", reduction.key = "UMAP")
Merged_Seurat_Rescaled <- FindNeighbors(Merged_Seurat_Rescaled, reduction = "harmony", dims= 1:30)
Merged_Seurat_Rescaled <- FindClusters(Merged_Seurat_Rescaled, resolution = 0.1, random.seed = SEED )

DimPlot(Merged_Seurat_Rescaled, label= TRUE)



#saveRDS(Merged_Seurat_Rescaled , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/3pops_Object.rds")







#If you want to use subclustering
Merged_Seurat_Rescaled2 = FindSubCluster(  Merged_Seurat_Rescaled,  cluster="4",
  graph.name= "RNA_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.2,
  algorithm = 1
)

Merged_Seurat_Rescaled2 = SetIdent(Merged_Seurat_Rescaled2, value= "sub.cluster")
DimPlot(Merged_Seurat_Rescaled2, label= TRUE)







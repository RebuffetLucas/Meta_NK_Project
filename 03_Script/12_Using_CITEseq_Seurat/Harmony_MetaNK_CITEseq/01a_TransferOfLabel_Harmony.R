## @knitr TransferLabelHarmony

#Make an object with  V2 , V3 and healthy  Tang PBMC
#Load all the objects
#Load V2 data
PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")
PBMC_Meta_NK_V2= SetIdent(PBMC_Meta_NK_V2, value= "FirstClust")

print(DimPlot(PBMC_Meta_NK_V2, cols = palette))


# Load Tang
PBMC_Blood_All = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Review_Nature/05_Output/CITEseq_Analysis/Seurat_3Clusters.rds")

#Pretreatment of the CITEseq data
PBMC_Blood_All  = subset(PBMC_Blood_All , subset = time == "0")

#PBMC_Blood_All@assays[["ADT"]] = NULL




print(DimPlot(PBMC_Blood_All, reduction = "wnn.umap"))




# table(PBMC_Blood_All$meta_histology)
# table(PBMC_Blood_All$datasets)
# 
# 
# table(droplevels(PBMC_Blood_All$meta_tissue))
# table(droplevels(PBMC_Blood_All$meta_histology))
# table(droplevels(PBMC_Blood_All$datasets))
# table(droplevels(PBMC_Blood_subset$datasets))
# table(droplevels(PBMC_Blood_subset$batch ))
# table(droplevels(PBMC_Blood_subset$sampleID))

#length(table(droplevels(PBMC_Blood_subset$sampleID)))



#Pick a subset of meta_histo to integrate (or not)
if(SUBSET_META_HISTO==TRUE){
PBMC_Blood_subset = subset(PBMC_Blood_All, subset = meta_histology ==  META_SUBSET)
}else{
PBMC_Blood_subset = PBMC_Blood_All
}


#Downsample randomly (or not)
if(SUBSET_DOWNSAMPLE==TRUE){
PBMC_Blood_subset = subset(PBMC_Blood_subset, downsample = N_DOWNSAMPLE)
}



#table(PBMC_Blood_subset$meta_histology)



#Keep genes intersection for each Dataset
genes_to_keep = intersect (rownames(PBMC_Meta_NK_V2)  , rownames(PBMC_Blood_subset ))

#Merge data
data_all =  merge(PBMC_Blood_subset, y = PBMC_Meta_NK_V2, add.cell.ids = c("Tang", "Meta_V2"))
data_all = subset(data_all, features= genes_to_keep)

data_all

#Meta Data
# metadata metaNK_V2 ####
metadata_NK_V2=PBMC_Meta_NK_V2@meta.data[,c("Dataset","orig.ident","FirstClust","SecondClust", "nkg2c")]
metadata_NK_V2$Project="MetaNK_V2"
metadata_NK_V2$meta_histology = as.factor("Healthy donor")
rownames(metadata_NK_V2) = paste0("Meta_V2_",rownames(metadata_NK_V2))
metadata_NK_V2=metadata_NK_V2[,
                              c("Dataset","orig.ident", "nkg2c", "Project","FirstClust","SecondClust", "meta_histology")]
metadata_NK_V2$Batch_Cor = metadata_NK_V2$orig.ident



# metadata Tang ####
metadata_Tang=PBMC_Blood_subset@meta.data[,c("meta_patientID","meta_histology","celltype","Majortype", "datasets", VARIABLE_TANG_FOR_BATCHCOR)]
colnames(metadata_Tang)[1]="orig.ident"
colnames(metadata_Tang)[6]="Batch_Cor"
metadata_Tang$FirstClust=metadata_Tang$Majortype
metadata_Tang$SecondClust=metadata_Tang$celltype
metadata_Tang$Dataset=metadata_Tang$datasets
metadata_Tang$Project="Tang"
metadata_Tang$nkg2c="unknown"
rownames(metadata_Tang) = paste0("Tang_",rownames(metadata_Tang))
metadata_Tang=metadata_Tang[,
                          c("Dataset","orig.ident", "nkg2c", "Project","FirstClust","SecondClust", "meta_histology", "Batch_Cor")]



#table(droplevels(metadata_Tang$Dataset))
#
#metadata_NK_V2[1:6, ]
#metadata_Tang[1:6, ]



# Combining metadata ####
metadata_all=rbind(metadata_NK_V2, metadata_Tang)
metadata_all=metadata_all[colnames(data_all),]

#Inserting them in the combined seurat object
data_all@meta.data = metadata_all


#table(data_all$FirstClust , data_all$Project)
#table(data_all$Dataset )
#table(droplevels( data_all$Batch_Cor ) )


#Normalize
data_all  = NormalizeData(data_all)


#Forcing variable Features:

if (VARIABLE_FEATURES_METHOD=="FORCE_VF1"){
  VF1 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VariableFeatures.rds")
  intersect(VF1, genes_to_keep)
  intersect(VariableFeatures(data_all) , intersect(VF1, genes_to_keep) )
  VARIABLE_GENES_INTERSECT = intersect(VF1, genes_to_keep)
  VariableFeatures(data_all)= VARIABLE_GENES_INTERSECT
  }

#Looking directly at the variable features of the global object:

if (VARIABLE_FEATURES_METHOD=="VARIABLE_FEATURE_BASIC"){
data_all = FindVariableFeatures(data_all, selection.method = "vst", nfeatures = 2000)
}

table( data_all$Dataset)
#Find Variable Features across projects
#Find variable features with Integration features across PROJECTS 
#Variable on orig.ident
if (VARIABLE_FEATURES_METHOD=="VARIABLE_FEATURE_INTEGRATION_FEATURES"){
seurat_list_Rescaled = SplitObject(data_all, split.by = "Project")
seurat_list_Rescaled <- lapply(seurat_list_Rescaled, function(x) FindVariableFeatures(x, nfeatures= 5000, selection.method = "vst"))
VARIABLE_FEATURES = SelectIntegrationFeatures(  seurat_list_Rescaled,  nfeatures = 2000,  assay = NULL,  verbose = TRUE)
}


#Scale and Run PCA
data_all = ScaleData(data_all, do.scale = DATA_SCALE, do.center = DATA_CENTER, features=  rownames(data_all) )
data_all= RunPCA(data_all, npcs = FINDCLUSTERS_USE_PCA_NBDIMS)



p1 = DimPlot(data_all, reduction= "pca", group.by = "Project")

#table(PBMC_Blood_All$datasets, PBMC_Blood_All$meta_histology)
#table(data_all$Project )
#table(droplevels( data_all$orig.ident ))

data_all$orig.ident = droplevels(data_all$orig.ident)
data_all$Batch_Cor = droplevels(data_all$Batch_Cor)

#table(data_all$Batch_Cor )



#Run Batch-Correction
if (HARMO_USE_REF ==TRUE){
data_all= data_all %>%                  
  RunHarmony("Batch_Cor", reference_values = HARMO_REF_TO_USE ,  plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:FINDCLUSTERS_USE_PCA_NBDIMS ) 
}


if (HARMO_USE_REF ==FALSE){
  data_all= data_all %>%                  
    RunHarmony("Batch_Cor",  plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:FINDCLUSTERS_USE_PCA_NBDIMS ) 
}



#Run UMAP
data_all <- RunUMAP(data_all, reduction = "harmony", dims= 1:FINDCLUSTERS_USE_PCA_NBDIMS, reduction.name = "UMAP", reduction.key = "UMAP")


#Clustering
data_all <- FindNeighbors(data_all, reduction = "harmony", dims= 1:FINDCLUSTERS_USE_PCA_NBDIMS)
data_all <- FindClusters(data_all, resolution = FINDCLUSTERS_RESOLUTION, random.seed = SEED )


#SaveRDS (or not)
if( SAVE_RDS==TRUE){
  saveRDS(data_all , RDS_DIR_OUTPUT)
}

#data_all$Dataset




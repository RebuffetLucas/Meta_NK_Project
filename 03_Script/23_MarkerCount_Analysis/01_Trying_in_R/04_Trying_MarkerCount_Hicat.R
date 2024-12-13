#Trying MarkerCount HiCAT version

#Silent starting messages

#suppressPackageStartupMessages(library(anndata))
suppressPackageStartupMessages(library(reticulate))

library(reticulate)

#load MarkerCount
mkrcnt <-  import("MarkerCount.hicat")


#Load data

## @knitr TransferOfLabel_StandardSeurat

#Make an object with  V2 , V3 and healthy  Tang PBMC
#Load all the objects
#Load V2 data
PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

#print(DimPlot(PBMC_Meta_NK_V2, cols = palette))


# Load Tang
PBMC_Blood_All = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5seurat")

#print(DimPlot(PBMC_Blood_All))



# table(PBMC_Blood_All$meta_histology)
# table(PBMC_Blood_All$datasets)
# 
# 
#table(droplevels(PBMC_Blood_All$meta_tissue))
# table(droplevels(PBMC_Blood_All$meta_histology))
# table(droplevels(PBMC_Blood_All$datasets))
# table(droplevels(PBMC_Blood_subset$datasets))
# table(droplevels(PBMC_Blood_subset$batch ))
# table(droplevels(PBMC_Blood_subset$sampleID))
# table(PBMC_Blood_All$meta_histology, PBMC_Blood_subset$datasets)
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

#data_all

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
metadata_all$Dataset= as.factor(metadata_all$Dataset)

#Inserting them in the combined seurat object
data_all@meta.data = metadata_all

table(data_all$Batch_Cor, data_all$Project)

table(data_all$Dataset )



#table(data_all$FirstClust , data_all$Project)
#table(data_all$Dataset )
#table(droplevels( data_all$Batch_Cor ) )


#Normalize
data_all_list = SplitObject(data_all, split.by = DATA_SPLIT) #May be we can change to not dataset but sample for building the ref

for (i in 1:length(data_all_list)) {
  data_all_list[[i]] <- NormalizeData(data_all_list[[i]], verbose = FALSE) #Maybe we need to change to MultiBatchNorm Here
  data_all_list[[i]] <- FindVariableFeatures(data_all_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

#Integration
reference.list = data_all_list[REFERENCE_FOR_INTEGRATION]
Integration.anchors= FindIntegrationAnchors(object.list = reference.list, dims= 1:30)

data_integrated =  IntegrateData(anchorset = Integration.anchors, dims = 1:30)

#Dataviz of the integrated
DefaultAssay(data_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
data_integrated <- ScaleData(data_integrated, verbose = FALSE)
data_integrated <- RunPCA(data_integrated, npcs = 30, verbose = FALSE)
data_integrated <- RunUMAP(data_integrated, reduction = "pca", dims = 1:30)



QUERY_DATASETS = setdiff(names(data_all_list) , REFERENCE_FOR_INTEGRATION)


#data.query = data_all_list[["GSE162025"]] #A transformer en une boucle en parourant tout les QUERY_DATASETS


#data_all_list[QUERY_DATASETS]
#data_all_list[QUERY_DATASETS][-1]

data.query = merge(data_all_list[QUERY_DATASETS][[1]] , y = data_all_list[QUERY_DATASETS][-1]) #Merge all the query together

table(data.query$meta_histology)


##############################################################################
######################### Using MarkerCount in R #############################
##############################################################################


#Load Data

adata_test <- data.query
adata_ref <- data_integrated

#Set data to pass to MarkerCount
X_ref = t(adata_ref@assays[["RNA"]]@counts)  # rows: cell, cols: genes
cell_type_ref = adata_ref@meta.data[["FirstClust"]]

X_test = t(adata_test@assays[["RNA"]]@counts)  # rows: cell, cols: genes

cell_type_test = adata_test@meta.data[["SecondClust"]]


#Check the levels
as.character(unique(cell_type_ref))
as.character(unique(cell_type_test))


#Run MarkerCount-Ref
df_res <- mkrcnt$MarkerCount_Ref( X_ref, cell_type_ref, 
                                  X_test = X_test, df_mkr_mat = NULL,
                                  N_mkrs = 18, of_th = 0.9, # min_th = 0.2, 
                                  #cell_types_to_excl = undesired, 
                                  cluster_label = NULL, X_pca = NULL, 
                                  log_transformed = FALSE, 
                                  file_to_save_marker = 'UC_cell_markers', 
                                  verbose = TRUE)

df_res[,'cell_type_org'] <- cell_type_test
head(df_res)




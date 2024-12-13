## @knitr TransferOfLabel_StandardSeurat

#Def function
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))}
  
  
  

#Make an object with  V2 , V3 and healthy  Tang PBMC
#Load all the objects
#Load V2 data
PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

#print(DimPlot(PBMC_Meta_NK_V2, cols = palette))

PBMC_Meta_NK_V2$SecondClust = PBMC_Meta_NK_V2$FirstClust
levels(PBMC_Meta_NK_V2$FirstClust)= c("NK1", "NK1", "NK1", "NK1", "NK2", "NK3")
PBMC_Meta_NK_V2 = SetIdent(PBMC_Meta_NK_V2, value = "FirstClust")
#DimPlot(PBMC_Meta_NK_V2, cols = palette2)

# Load Tang
PBMC_Blood_All = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5seurat")
#PBMC_Blood_All = LoadH5Seurat("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5seurat")

#table(PBMC_Blood_All$meta_tissue , PBMC_Blood_All$meta_histology)
#print(DimPlot(PBMC_Blood_All))


#table(PBMC_Blood_All$datasets, PBMC_Blood_All$meta_histology)
#table(PBMC_Blood_All$meta_histology)


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


if(SUBSET_META_TISSUE==TRUE){
  PBMC_Blood_subset = subset(PBMC_Blood_All, subset = meta_tissue %in%  META_TISSUE_SUBSETS)
  PBMC_Blood_subset$meta_tissue = droplevels(PBMC_Blood_subset$meta_tissue )
}else{
  PBMC_Blood_subset = PBMC_Blood_All
}

#Pick a subset of meta_histo to integrate (or not)
if(SUBSET_META_HISTO==TRUE){
PBMC_Blood_subset = subset(PBMC_Blood_All, subset = meta_histology ==  META_SUBSET)
}else{
PBMC_Blood_subset = PBMC_Blood_subset
}

#Downsample randomly (or not)
if(SUBSET_DOWNSAMPLE==TRUE){
PBMC_Blood_subset = subset(PBMC_Blood_subset, downsample = N_DOWNSAMPLE)
}



#Create the tables for after filtering
table_samples = data.frame( table(PBMC_Blood_subset$meta_tissue , PBMC_Blood_subset$meta_histology))
colnames(table_samples) = c("Tissue" , "Cancer Type", "Number of Cells")
table_All2 = datatable(table_samples)

table_samples_Light = table_samples[table_samples$`Number of Cells`!=0,]
table_Light2 = datatable(table_samples_Light)



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
metadata_NK_V2$meta_tissue= as.factor("blood")
metadata_NK_V2$Batch_Cor = metadata_NK_V2$orig.ident




# metadata Tang ####
metadata_Tang=PBMC_Blood_subset@meta.data[,c("meta_patientID","meta_histology","celltype","Majortype", "meta_tissue","datasets", DATA_SPLIT)]
colnames(metadata_Tang)[1]="orig.ident"
colnames(metadata_Tang)[length(colnames(metadata_Tang))] = "Batch_Cor"
metadata_Tang$FirstClust=metadata_Tang$Majortype
metadata_Tang$SecondClust=metadata_Tang$celltype
metadata_Tang$Dataset=metadata_Tang$datasets
metadata_Tang$Project="Tang"
metadata_Tang$nkg2c="unknown"
rownames(metadata_Tang) = paste0("Tang_",rownames(metadata_Tang))
metadata_Tang=metadata_Tang[,
                          c("Dataset","orig.ident", "nkg2c", "Project","FirstClust","SecondClust", "meta_histology", "meta_tissue", "Batch_Cor")]



if(DO_USE_COMPOSED_DATA_SPLIT_FOR_TANG ==TRUE){
  metadata_Tang[,"Batch_Cor"] = paste0(metadata_Tang[,VAR_SPLIT1] ,"_", metadata_Tang[,VAR_SPLIT2])
}


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

#table(data_all$Batch_Cor, data_all$Project)
#table(data_all$Dataset )



#table(data_all$FirstClust , data_all$Project)
#table(data_all$Dataset )
#table(droplevels( data_all$Batch_Cor ) )




if( DO_USE_CATEGORY_TRESHOLD ==TRUE){
  Batch_Cor_keep = names(table(data_all$Batch_Cor)[table(data_all$Batch_Cor) > MINIMAL_NUMBER_BATCH_COR_CATEGORY])
  data_all2 = subset(data_all, subset = Batch_Cor %in% Batch_Cor_keep  )
}



data_all_list = SplitObject(data_all2, split.by = "Batch_Cor") #May be we can change to not dataset but sample for building the ref


#Normalize
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

data_integrated$FirstClust = as.factor(data_integrated$FirstClust )

#data.query = data_all_list[["GSE162025"]] #A transformer en une boucle en parourant tout les QUERY_DATASETS


#data_all_list[QUERY_DATASETS]
#data_all_list[QUERY_DATASETS][-1]


data.query = merge(data_all_list[QUERY_DATASETS][[1]] , y = data_all_list[QUERY_DATASETS][-1]) #Merge all the query together


#Projecting on the integrated reference

data.anchors <- FindTransferAnchors(reference = data_integrated, query = data.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = data.anchors, refdata = data_integrated$FirstClust , 
                            dims = 1:30)

data.query <- AddMetaData(data.query, metadata = predictions)

data.query$tissue_Firstclust = paste0(data.query$meta_tissue , "_", data.query$SecondClust )
data.query$tissue_Histo = paste0(data.query$meta_tissue , "_", data.query$meta_histology )

#table(data.query$predicted.id )

if(DO_SAVE_QUERYandREF ==TRUE){
  
saveRDS(data.query, DATA_QUERY_SAVE_PATH)
saveRDS(data_integrated, DATA_REF_SAVE_PATH)

}

## @knitr TransferOfLabel_StandardSeurat

#Def function
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))}
  
  
  

#Make an object with  V2 , V3 and healthy  Tang PBMC
#Load all the objects
#Load V2 data
PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

PBMC_Meta_NK_V2$SecondClust = PBMC_Meta_NK_V2$FirstClust
levels(PBMC_Meta_NK_V2$FirstClust)= c("NK1", "NK1", "NK1", "NK1", "NK2", "NK3")

#print(DimPlot(PBMC_Meta_NK_V2, cols = palette))
#table(PBMC_Meta_NK_V2$FirstClust)

df = PBMC_Meta_NK_V2@meta.data 
df$row_names <- rownames(df)  # If your data doesn't have row names, first assign them

set.seed(SEED)

# Downsample 20% of the rows for each FirstClust category
downsampled_df <- ddply(df, .(FirstClust), function(sub_df) {
  # Calculate number of rows to sample
  sample_size <- ceiling(nrow(sub_df) * PERCENT_DOWNSAMPLE)
  
  # Sample row names
  sampled_rows <- sample(sub_df$row_names, size = sample_size)
  
  # Return downsampled subset
  sub_df[sub_df$row_names %in% sampled_rows, ]
})

#table(downsampled_df$FirstClust)

#downsampled_df$row_names

PBMC_Blood_subset = subset(PBMC_Meta_NK_V2, cells = downsampled_df$row_names)


#table(PBMC_Blood_subset$orig.ident)

#DimPlot(PBMC_Meta_NK_V2, cols = palette2)

# Load Tang
#PBMC_Blood_All = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5seurat")



PBMC_Meta_NK_V2 = subset(PBMC_Meta_NK_V2, cells= Cells(PBMC_Blood_subset), invert = TRUE)
PBMC_Meta_NK_V2 = SetIdent(PBMC_Meta_NK_V2, value = "FirstClust")

PBMC_Blood_subset$Projet = "Query"
PBMC_Meta_NK_V2$Projet = "Ref"
#Keep genes intersection for each Dataset
genes_to_keep = intersect (rownames(PBMC_Meta_NK_V2)  , rownames(PBMC_Blood_subset ))

#Merge data
data_all =  merge(PBMC_Blood_subset, y = PBMC_Meta_NK_V2, add.cell.ids = c("Query", "Ref"))
data_all = subset(data_all, features= genes_to_keep)


#table(data_all$Batch_Cor, data_all$Project)
#table(data_all$Dataset )

data_all$split = paste0(data_all$orig.ident ,  data_all$Projet)

#table(data_all$FirstClust , data_all$Project)
#table(data_all$Dataset )
#table(droplevels( data_all$Batch_Cor ) )

#Normalize
data_all_list = SplitObject(data_all, split.by = "orig.ident") #May be we can change to not dataset but sample for building the ref

for (i in 1:length(data_all_list)) {
  data_all_list[[i]] <- NormalizeData(data_all_list[[i]], verbose = FALSE) #Maybe we need to change to MultiBatchNorm Here
  data_all_list[[i]] <- FindVariableFeatures(data_all_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


#Re-split by query and orig ident
data_all2 = merge(data_all_list[[1]] , y = data_all_list[-1]) #Merge all the query together

data_all_list = SplitObject(data_all2, split.by = "split") #May be we can change to not dataset but sample for building the ref



#Integration
REFERENCE_FOR_INTEGRATION = names(data_all_list)[grepl("Ref" , names(data_all_list) )]
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

data.query <- AddMetaData(data.query, metadata = predictions) #Predictions 4 pops


predictions2 <- TransferData(anchorset = data.anchors, refdata = data_integrated$SecondClust , 
                            dims = 1:30)

data.query2 <- AddMetaData(data.query, metadata = predictions2) #Predictions 6 pops

#table(data.query$predicted.id )


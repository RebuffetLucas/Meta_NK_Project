
METADATA_QUERY = "wsnn_res.0.5"
remotes::install_github("carmonalab/SignatuR")
library(SignatuR)

library(scuttle)

#Upload differencially expressed Genes
DEG_DESQ2 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/230626CVa_markers_vf1/resFirstClust.rds")
DEG_Seurat = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/230626CVa_markers_vf1/resFirstClust_Seurat.rds")


PBMC = Tang_Object2
PBMC=NormalizeData(PBMC)

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

Score_Level1 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds")
Score_6Clusters = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP11_6_Clusters.rds")
Score_5Clusters = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Merging_NK1AandNK1B/TOPAllupto68_5MainClusters.rds")


Score_Malmberg= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg.rds")



MONITORED_Markers  = Score_5Clusters
#Pre-handling ONLY FOR MALMBERGDATA !!!!
MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(PBMC)))
min_genes= min(unlist(lapply(MONITORED_Markers, FUN=function(x) length(x))))
MONITORED_Markers = lapply(MONITORED_Markers, FUN = function(x) head(x, n=min_genes) )


#Go with Direct scoring
for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = MONITORED_Markers[[i]], pool= NULL ,name= i , seed=19)
}




#Remove the "1" at the End of the modulescore
#colnames(PBMC@meta.data) <-
#  gsub(x = colnames(PBMC@meta.data)
#       , pattern = "1$"
#       , replacement = ""
#  )

#Print a Pheatmap of the scoring at the level of the clusters ( METADATA_QUERY)
table = PBMC@meta.data[,c( as.character(METADATA_QUERY) , paste0(names(MONITORED_Markers),1) )]
table$group_of_interest = table[,as.character(METADATA_QUERY)]
table[,as.character(METADATA_QUERY)]=NULL

df = table %>% 
  group_by(group_of_interest) %>% 
  dplyr::summarise(across(everything(), list(mean)))

df= as.data.frame(df)
rownames(df)=df$group_of_interest

df2= df[,-1]

    #Print
pheatmap(as.matrix(df2), scale= "column", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_col = 15, fontsize_row = 15)
pheatmap(as.matrix(df2), scale= "column")


  #Direct assignation on the basis of AddModuleScore
table2 = table[, colnames(table) != 'group_of_interest'] %>% mutate(Class=names(.)[max.col(.)])
table2$group_of_interest = table$group_of_interest
table(table2$Class, table2$group_of_interest)

#Look at scoring on a FeaturePlot
FeaturePlot(PBMC, features= paste0(names(MONITORED_Markers),1) , reduction="wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



 #Trying assignement using HTO Does not work because  it needs an 

  
  #Build the Data slot based on ModuleScores
table_data_HTO = as.matrix(table[, colnames(table) != 'group_of_interest'])
table_data_HTO= t(table_data_HTO)
table_data_HTO[,1:5]
table_count_HTO[,1:5]


  #Built the count slot based on counts of genes composing the modulescore
table_count_HTO = table_data_HTO

#i = names(MONITORED_Markers)[1]
#sum(rownames(PBMC_counts) %in% as.character(unlist(MONITORED_Markers[[i]])))

for (i in names(MONITORED_Markers)){
print(i)
ids = rownames(PBMC_counts)
ids[ids%in% as.character(unlist(MONITORED_Markers[[i]]))]= paste0(i, 1)
ids[ids!= paste0(i, 1)]= NA
table_count_HTO2 = sumCountsAcrossFeatures(x= PBMC@assays[["RNA"]]@counts,ids = ids )
print(table_count_HTO2[,1:5])
table_count_HTO[paste0(i,1),]=table_count_HTO2
print(table_count_HTO[paste0(i,1),1:5])
}



pbmc.hashtag= PBMC
pbmc.hashtag= FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

  # Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(data = table_data_HTO)
pbmc.hashtag[["HTO"]]@counts=table_count_HTO

dim(pbmc.hashtag[["HTO"]]@counts)
dim(pbmc.hashtag[["HTO"]]@data)

pbmc.hashtag[["HTO"]]@counts[,1:10]
pbmc.hashtag[["HTO"]]@data[,1:10]


#Demultiplexing
      pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

#Visualize demultiplexing results

table(pbmc.hashtag$HTO_maxID, PBMC$FirstClust)
table(pbmc.hashtag$HTO_maxID)
table(pbmc.hashtag$HTO_secondID)
table(pbmc.hashtag$HTO_classification.global)
levels(pbmc.hashtag$HTO_classification.global)

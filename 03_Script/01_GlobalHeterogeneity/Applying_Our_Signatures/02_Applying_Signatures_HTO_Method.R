
METADATA_QUERY = "celltype.l3"

library(scuttle)

#Upload differencially expressed Genes
DEG_DESQ2 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/230626CVa_markers_vf1/resFirstClust.rds")
DEG_Seurat = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/230626CVa_markers_vf1/resFirstClust_Seurat.rds")


PBMC = Tang_Object2
PBMC=NormalizeData(PBMC)

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")


#Extract Top20
DEG_Seurat %>% 
  group_by(cluster) %>%
  filter(gene %in% rownames(PBMC) ) %>%
  filter(p_val_adj<0.05 & avg_log2FC>0 ) %>%
  top_n(n=20 , wt = avg_log2FC) -> top10


table(top10$cluster)


List_Top10=list()

for (i in levels(top10$cluster)){
  top10 %>% 
    filter(cluster== i ) -> top10_cluster
  print(top10_cluster$gene)
  List_Top10  = append(List_Top10, data.frame(top10_cluster$gene))
}

names(List_Top10) = levels(top10$cluster)

intersect(List_Top10$NK2A, List_Top10$NK2B)

List_Top10= lapply(List_Top10, as.data.frame)

#Change colnames
for (i in names(List_Top10)){
  print(i)
  colnames(List_Top10[[i]]) = i
}

MONITORED_Markers  = List_Top10

#DotPlot(PBMC, features = unique(unlist(List_Top10)) , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))


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

#Vizualisation at sc level
PBMC$max_identity = table2$Class
DimPlot(PBMC, group.by = "max_identity", pt.size= 0.4)
DimPlot(PBMC)


FeaturePlot(PBMC, features = levels(as.factor(table2$Class)), pt.size= 0.6)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Bar graph
p5 = ggplot(PBMC@meta.data, aes(x=FirstClust, fill= max_identity)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(PBMC@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p7 = ggplot(PBMC@meta.data, aes(x=Dataset, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  

p5


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
ids = rownames(PBMC@assays[["RNA"]]@counts)
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
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR") # Useless step (just to make sure that Normalized data becomes the active slot)
pbmc.hashtag[["HTO"]]@data = table_data_HTO

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

#Go Dataviz with RidgePlot

Idents(pbmc.hashtag) = "HTO_maxID"

RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:6], ncol = 2, group.by = "FirstClust" )

p1 = FeatureScatter(pbmc.hashtag, feature1 = "NK1A1", feature2 = "NK1B1", group.by = "FirstClust")
p2 = FeatureScatter(pbmc.hashtag, feature1 = "NK2A1", feature2 = "NK2B1", group.by = "FirstClust")
p3 = FeatureScatter(pbmc.hashtag, feature1 = "NK31", feature2 = "NK1C1", group.by = "FirstClust")
p4 = FeatureScatter(pbmc.hashtag, feature1 = "NK2B1", feature2 = "NK1A1", group.by = "FirstClust")

p1+p2+p3+p4


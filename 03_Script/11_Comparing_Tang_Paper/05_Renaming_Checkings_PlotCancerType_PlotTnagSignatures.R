#Implementing our nomenclature / checking the signatures / Producing barplot across cancer type / comparing with their nomenclature

data_all = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/Integrated_Seurat_Object/IntegratedMetaV2_15Clusters_Tang.rds")

DimPlot(data_all)

#Get a pie chart of the proportions of each cluster in the pop

df= as.data.frame(table(Idents(data_all)))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

p6 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    
p6 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    



PBMC= subset( data_all, idents= c("15", "14"), invert= TRUE )

p1 = DimPlot(PBMC, group.by = "FirstClust", split.by = "Project")
p2 = DimPlot(PBMC, group.by = "SecondClust", split.by = "Project")

p1 / p2

PBMC$Rename_Group = PBMC$RNA_snn_res.0.8


#renaming
PBMC$Rename_Group = droplevels(PBMC$Rename_Group)
levels(PBMC$Rename_Group) = c("NK1A", "NK1B", "NK1B", "NK3", "NK2B", "NK3", "NK1B", "NK3", "NK1C", "NK3", "NK1C", "NK2B", "NK2A", "NK3")

PBMC$Rename_Group = factor(PBMC$Rename_Group , levels = c("NK1C" , "NK3","NK1A", "NK1B",   "NK2B", "NK2A"))

levels(PBMC$Rename_Group) = relevel(levels(PBMC$Rename_Group), "NK3")


PBMC$Rename_Main_Group = droplevels(PBMC$Rename_Group)
levels(PBMC$Rename_Main_Group) = c("NK1", "NK1", "NK3", "NK2", "NK1", "NK2")



#Plots
DimPlot(PBMC, group.by = "Rename_Group", split.by = "Project")
DimPlot(PBMC, group.by = "Rename_Main_Group", split.by = "Project")

#Proportions Plot

PBMC$meta_histology= droplevels(PBMC$meta_histology)


#Bar Diagram of Datasets, Chemistry , Orig.ident:

PBMC$seurat_clusters = PBMC$Rename_Group
p5 = ggplot(PBMC@meta.data, aes(x=meta_histology, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=15),axis.text.x = element_text(angle = 70, hjust = 1)  ) 

p5

PBMC$seurat_clusters = PBMC$Rename_Main_Group
p5 = ggplot(PBMC@meta.data, aes(x=meta_histology, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90), axis.text.align = 1)  
p5 = ggplot(PBMC@meta.data, aes(x=meta_histology, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=15),axis.text.x = element_text(angle = 70, hjust = 1)  ) 

p5

table(PBMC$meta_histology)


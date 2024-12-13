#Plot and diagnostic of a clustering

BiocManager::install("limma")
#Read Seurat Object


Merged_Seurat_Rescaled = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK4/05_Output/01_GlobalHeterogeneity/PBMC_All.rds")


#BiocManager::install("DESeq2")

library(DESeq2)


#BiocManager::install("limma")
library(limma)

  
#Global vizualization of the UMAP
p1 = DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", label= TRUE, label.size = 5)
p2= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "Dataset")
p3= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "orig.ident")
p4= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "Chemistry")

p1
p1 + p2 + p3 + p4

p5= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "nkg2c")

DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", label = TRUE, label.size = 7)

#Get a pie chart of the proportions of each cluster in the pop

df= as.data.frame(table(Idents(Merged_Seurat_Rescaled)))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

p6 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    
p6 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    

p5 + p6 

p5= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "Firstclust")
p12= VlnPlot(Merged_Seurat_Rescaled, feature = "percent.mito")



p10= FeaturePlot(Merged_Seurat_Rescaled, features = c("FGFBP2", "FCGR3A", "SPON2", "GZMB"))
p11 = FeaturePlot(Merged_Seurat_Rescaled, features = c("XCL1", "XCL2", "GZMK", "IL7R"))

p11 = FeaturePlot(Merged_Seurat_Rescaled, features = c("CCL4", "CCL3"))


p14= FeaturePlot(Merged_Seurat_Rescaled, features = c("CD52" , "CCL5", "KLRC2", "IL32"))
p15 = FeaturePlot(Merged_Seurat_Rescaled, features = "PTGDS")


p12= FeaturePlot(Merged_Seurat_Rescaled, feature = "percent.mito", pt.size = 0.6)
p13= FeaturePlot(Merged_Seurat_Rescaled, feature = "percent.ribo")

p12= VlnPlot(Merged_Seurat_Rescaled, feature = "percent.mito", pt.size = 0)
p13= VlnPlot(Merged_Seurat_Rescaled, feature = "percent.ribo", pt.size = 0)


VlnPlot(Merged_Seurat_Rescaled, feature = "CXCR4", group.by = "Dataset")


p10
p11
p14
p15
p12 +p13

#Bar Diagram of Datasets, Chemistry , Orig.ident:

Merged_Seurat_Rescaled$seurat_clusters = Merged_Seurat_Rescaled$FirstClust
p5 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=orig.ident, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p7 = ggplot(Merged_Seurat_Rescaled@meta.data, aes(x=Dataset, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  

p5 + p6  / p7

#Global Diagnostic:
  #Recalculate
nCount = colSums(x = Merged_Seurat_Rescaled, slot = "counts")  # nCount_RNA
nFeature = colSums(x = GetAssayData(object = object, slot = "counts") > 0)  # nFeatureRNA

  #Plot
p9 = VlnPlot(Merged_Seurat_Rescaled, features = c( "nCount_RNA" , "nFeature_RNA"  ,"percent.ribo" ,"percent.mito"), group.by = "orig.ident", pt.size = 0)
p10 = VlnPlot(Merged_Seurat_Rescaled, features = c( "nCount_RNA" , "nFeature_RNA"  ,"percent.ribo" ,"percent.mito"),  pt.size = 0)

DotPlot(Merged_Seurat_Rescaled, features = c( "BAG3" , "NR4A1"  ,"FOS" ,"JUN"))
DotPlot(Merged_Seurat_Rescaled, features = c( "NEAT1"))


p9
p10


#Look at markers previously calculated
#Find Markers
All_Markers = FindAllMarkers(Merged_Seurat_Rescaled , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


#Only for a pop
Markers1 = FindMarkers(Merged_Seurat_Rescaled , ident.1 = "2", ident.2 = "4", only.pos = TRUE)

#Look at markers like Constance
markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Merged_Seurat_Rescaled)))) {
  print(cl)
  markers_overCl_cl=markers_overCl[
    markers_overCl$avg_log2FC>0&
      markers_overCl$p_val_adj<=0.05&
      markers_overCl$cluster==cl,
  ]
  
  markers_overCl_cl=	
    markers_overCl_cl[order(markers_overCl_cl$p_val),]
  topGenes_cl=markers_overCl_cl$gene[
    1:min(FINDMARKERS_SHOWTOP,nrow(markers_overCl_cl))]
  topGenes[[cl]]=topGenes_cl
  print(	topGenes_cl)
}

topGenes2

saveRDS(topGenes2, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/signatures/top20.rds")

top = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/signatures/top20.rds")


topGenesb=unique(unlist(topGenes))

DotPlot(Merged_Seurat_Rescaled, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

DotPlot(Merged_Seurat_Rescaled, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

DotPlot(Merged_Seurat_Rescaled, features = c("SPON2", "GZMB", "FCGR3A", "FGFBP2") , cols = "Spectral") + NoLegend() + theme(axis.text.x = element_text(angle = 90))












#Look the marker of a specific subset
All_Markers %>%
  filter(cluster== "GammaDelta") %>%
  slice_max(n=100, order_by = avg_log2FC)


All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) -> top10


DotPlot(Merged_Seurat_Rescaled, features = unique(top10$gene) , cols = "Spectral") + NoLegend() + theme(axis.text.x = element_text(angle = 90))


All_Markers %>%
  group_by(cluster) %>%
  filter(p_val_adj<=0.05 & avg_log2FC>0) %>%
  slice(n = 1:15) -> top10





#Representing it with pheatmap
cluster.averages <- AverageExpression(Merged_Seurat_Rescaled, group.by = "Secondclust", features = unique(top10$gene), slot= "data")
#cluster.averages$Secondclust = names(cluster.averages@active.ident)
pheatmap(cluster.averages$RNA, scale= "row", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_col = 15, fontsize_row = 10)














#Find Markers with DESEQ2

All_Markers = FindAllMarkers(Merged_Seurat_Rescaled , test.use= "DESeq2" )

saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/01_GlobalHeterogeneity/Output_First_Clustering/ChemData100/ChemData_Rescaled_Norm100/All_Markers_NK1_DESeq2.rds")

Markers7 = FindMarkers(Merged_Seurat_Rescaled , ident.1 = "7")

#Look the marker of a specific subset
All_Markers %>%
  filter(cluster=="NK_2B") %>%
  slice_max(n=100, order_by = avg_log2FC)








#Look at their automatic assignement

sce <- as.SingleCellExperiment(Merged_Seurat_Rescaled)
BlueprintEncodeData.ref <- celldex::BlueprintEncodeData()
BlueprintEncodeData.fine <- SingleR(test = sce,ref = BlueprintEncodeData.ref,labels = BlueprintEncodeData.ref$label.fine)
Merged_Seurat_Rescaled@meta.data$BlueprintEncodeData.fine <- BlueprintEncodeData.fine$labels

DimPlot(Merged_Seurat_Rescaled,group.by = "BlueprintEncodeData.fine")




#Save Object
  #For NK_All
#saveRDS(All_Markers , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/01_GlobalHeterogeneity/Output_First_Clustering/ChemData100/ChemData_Rescaled_Norm100/NK_All_Top_Markers_15clusters.rds")
#saveRDS(Merged_Seurat_Rescaled , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/01_GlobalHeterogeneity/Output_First_Clustering/ChemData100/ChemData_Rescaled_Norm100/NK_All_15clusters.rds")

saveRDS( All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/01_GlobalHeterogeneity/Output_First_Clustering/ChemData100/ChemData_Rescaled_Norm100/send_to_constance/Markers_NK1_All.rds")






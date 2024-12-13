#Identify diffrerent subpops:

library(dplyr)

#For Dim
Merged_Seurat_Rescaled$FirstClust =as.factor(Merged_Seurat_Rescaled$FirstClust)
Merged_Seurat_Rescaled$SecondClust =as.character(Merged_Seurat_Rescaled$FirstClust)

df = as.data.frame(Merged_Seurat_Rescaled$FirstClust)
colnames(df) = "subcluster"

Merged_Seurat_Rescaled2=  FindSubCluster(
  Merged_Seurat_Rescaled,
  cluster= "0",
  graph.name= "RNA_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.2,
  algorithm = 1
)


Merged_Seurat_Rescaled2=SetIdent(Merged_Seurat_Rescaled2, value= "sub.cluster")
DimPlot(Merged_Seurat_Rescaled2, group.by = "sub.cluster")

#Retrieve Cell Names
List_mNK1 = Cells(subset(Merged_Seurat_Rescaled2, idents="0_0" ))
List_mNK2 = Cells(subset(Merged_Seurat_Rescaled2, idents="0_1" ))


#For Adapt
Merged_Seurat_Rescaled2=  FindSubCluster(
  Merged_Seurat_Rescaled,
  cluster= "1",
  graph.name= "RNA_snn",
  subcluster.name = "sub.cluster",
  resolution = 0,
  algorithm = 1
)

Merged_Seurat_Rescaled2=SetIdent(Merged_Seurat_Rescaled2, value= "sub.cluster")
DimPlot(Merged_Seurat_Rescaled2, group.by = "sub.cluster")

List_aNK1 = Cells(subset(Merged_Seurat_Rescaled2, idents="1_0" ))



#For Bright
Merged_Seurat_Rescaled2=  FindSubCluster(
  Merged_Seurat_Rescaled,
  cluster= "2",
  graph.name= "RNA_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.2,
  algorithm = 1
)

Merged_Seurat_Rescaled2=SetIdent(Merged_Seurat_Rescaled2, value= "sub.cluster")
DimPlot(Merged_Seurat_Rescaled2, group.by = "sub.cluster")

List_iNK1 = Cells(subset(Merged_Seurat_Rescaled2, idents="2_0" ))
List_iNK2 = Cells(subset(Merged_Seurat_Rescaled2, idents="2_1" ))


#Rename
Merged_Seurat_Rescaled$SecondClust[List_iNK1]= "iNK1"
Merged_Seurat_Rescaled$SecondClust[List_iNK2]= "iNK2"

Merged_Seurat_Rescaled$SecondClust[List_mNK1]= "mNK1"
Merged_Seurat_Rescaled$SecondClust[List_mNK2]= "mNK2"


Merged_Seurat_Rescaled$SecondClust[List_aNK1]= "aNK"

Merged_Seurat_Rescaled$SecondClust= as.factor(Merged_Seurat_Rescaled$SecondClust)
Merged_Seurat_Rescaled= SetIdent(Merged_Seurat_Rescaled, value= "SecondClust")

DimPlot(Merged_Seurat_Rescaled) + scale_color_manual(values=palette)


#Look at markers previously calculated
#Find Markers

All_Markers = FindAllMarkers(Merged_Seurat_Rescaled2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)





All_Markers = FindAllMarkers(Merged_Seurat_Rescaled , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


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

topGenes2=topGenes

topGenesb=unique(unlist(topGenes))

DotPlot(Merged_Seurat_Rescaled, features = topGenesb , cols = "Spectral") + NoLegend() + theme(axis.text.x = element_text(angle = 90))

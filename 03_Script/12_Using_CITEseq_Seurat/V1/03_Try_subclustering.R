#Try subclustering of the CITE-seq Dataset


Seurat_NK2 = subset(Seurat_NK, idents= "NK Proliferating" , invert= TRUE)

#Noramlize both ADT and RNA


#Multimodal_Clustering
Seurat_NK2 <- FindMultiModalNeighbors(  Seurat_NK2, reduction.list = list("pca", "apca"),   dims.list = list(1:50, 1:50), modality.weight.name = "RNA.weight" )


Seurat_NK2 <- RunUMAP(Seurat_NK2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_NK2 <- FindClusters(Seurat_NK2, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE)

p1 <- DimPlot(Seurat_NK2, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) #+ NoLegend()
p2 <- DimPlot(Seurat_NK2, reduction = 'wnn.umap', group.by = 'celltype.l3', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2


DimPlot(Seurat_NK2, reduction= "wnn.umap")


#Look at signatures
#have a look at their signatures
DefaultAssay(object = Seurat_NK2) <- "SCT"
All_Markers = FindAllMarkers(Seurat_NK2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level

markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Seurat_NK2)))) {
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

DotPlot(Seurat_NK2, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Seurat_NK2, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Look at markers at the protein level

DefaultAssay(object = Seurat_NK2) <- "ADT"

All_Markers = FindAllMarkers(Seurat_NK2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level

markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Seurat_NK2)))) {
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

DotPlot(Seurat_NK2, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Seurat_NK2, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))




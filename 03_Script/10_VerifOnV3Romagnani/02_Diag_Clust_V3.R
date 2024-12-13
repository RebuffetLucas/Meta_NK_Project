#Diagnostic Verif V3 Romagnani

p1 = DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", label= TRUE, label.size = 5)
p2= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "nkg2c")
p3= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "orig.ident")
p4= DimPlot(Merged_Seurat_Rescaled, reduction = "UMAP", group.by = "Chemistry")

p2 + p3
p1 + p2 + p3 + p4



p10= FeaturePlot(Merged_Seurat_Rescaled, features = c("FGFBP2", "FCGR3A", "SPON2", "GZMB"))
p11 = FeaturePlot(Merged_Seurat_Rescaled, features = c("XCL1", "XCL2", "GZMK", "IL7R"))
p12 = FeaturePlot(Merged_Seurat_Rescaled, features = c("CCL4", "CCL3"))
p13= FeaturePlot(Merged_Seurat_Rescaled, features = c("CD52" , "CCL5", "KLRC2", "IL32"))

p10
p11
p12
p13




#Look at markers previously calculated
#Find Markers
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

DotPlot(Merged_Seurat_Rescaled, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

DotPlot(Merged_Seurat_Rescaled, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

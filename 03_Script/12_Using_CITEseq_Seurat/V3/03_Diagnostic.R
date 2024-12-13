# Vizualization

p1 = DimPlot(Seurat_NK_Final, reduction= "wnn.umap")
p2 = DimPlot(Seurat_NK_Final, reduction= "wnn.umap" , group.by = "celltype.l3")
p3 = DimPlot(Seurat_NK_Final, reduction= "wnn.umap" , group.by = "orig.ident")

p1 + p2 + p3



#Have a look at their signatures
DefaultAssay(object = Seurat_NK_Final) <- "SCT"

All_Markers = FindAllMarkers(Seurat_NK_Final , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level


markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Seurat_NK_Final)))) {
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

DotPlot(Seurat_NK_Final, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Seurat_NK_Final, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Look at markers at the protein level

DefaultAssay(object = Seurat_NK_Final) <- "ADT"
All_Markers = FindAllMarkers(Seurat_NK_Final , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level

markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Seurat_NK_Final)))) {
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

DotPlot(Seurat_NK_Final, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Seurat_NK_Final, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Look at markers at the protein level

DefaultAssay(object = Seurat_NK_Final) <- "ADT"

FeaturePlot(Seurat_NK_Final, features = "CD16", max.cutoff =  4)  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(Seurat_NK_Final, features = "CD56-1", max.cutoff =  4)  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Diagnostic of the clustering


#Seurat_NK = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")

#Have a look at the signatures
#have a look at their signatures
All_Markers = FindAllMarkers(Seurat_NK , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level

DefaultAssay(object = Seurat_NK) <- "SCT"

markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Seurat_NK)))) {
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

DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Look at markers at the protein level

DefaultAssay(object = Seurat_NK) <- "ADT"
All_Markers = FindAllMarkers(Seurat_NK , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level

markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(Seurat_NK)))) {
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

DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Look at markers at the protein level

DefaultAssay(object = Seurat_NK) <- "ADT"


FeaturePlot(Seurat_NK, features= c("CD16"), reduction = "wnn.umap")  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

DimPlot(Seurat_NK, group.by = "orig.ident", reduction = "wnn.umap")  



Seurat_NK = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")

Seurat_NK= subset(Seurat_NK, idents= "NK Proliferating", invert=TRUE)

DefaultAssay(object = Seurat_NK) <- "ADT"
Seurat_NK = NormalizeData(Seurat_NK, normalization.method = 'CLR', margin = 2)

DimPlot(Seurat_NK, reduction = "wnn.umap")
#Have a look at the signatures
#have a look at their signatures


DefaultAssay(object = Seurat_NK) <- "SCT"
All_Markers = FindAllMarkers(Seurat_NK , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/HAO_Clusters/All_Markers_HAO_5clusters.rds")
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


p1 = DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

p1
#DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Save the figure
ggsave("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3Clusters_DotPlot_RNA.png",
       plot = p1 , width = 15, height = 35, dpi = 600, units = "cm")

#Save the differential expression


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

p1 = DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral")  + coord_flip() + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 8)) 
p1
#DotPlot(Seurat_NK, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

#Save the figure
ggsave("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3Clusters_DotPlot_ADT.png",
       plot = p1 , width = 15, height = 20, dpi = 600, units = "cm")


VlnPlot(Seurat_NK, features= "CD160", pt.size = 0)


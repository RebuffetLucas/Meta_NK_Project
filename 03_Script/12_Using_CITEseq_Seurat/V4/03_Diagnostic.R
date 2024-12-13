# Vizualization



Seurat_NK_Final =readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_After_Satija_Procedure_D0ref_D3_D7_Optimized.rds")
p1 = DimPlot(Seurat_NK_Final, reduction= "wnn.umap", label = TRUE, label.size = 6)
p2 = DimPlot(Seurat_NK_Final, reduction= "wnn.umap" , group.by = "celltype.l3", label = TRUE)
p3 = DimPlot(Seurat_NK_Final, reduction= "wnn.umap" , group.by = "orig.ident", label = F)

p1 + p2 + p3

#BarGrpah for Diag

df= as.data.frame(table(Idents(Seurat_NK_Final)))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

p4 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    


p5 = ggplot(Seurat_NK_Final@meta.data, aes(x=orig.ident, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Seurat_NK_Final@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p7 = ggplot(Seurat_NK_Final@meta.data, aes(x=donor, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  

p4 + p5 + p6  + p7

Seurat_NK_Final[["percent.mt"]] <- PercentageFeatureSet(Seurat_NK_Final, pattern = "^MT-")
VlnPlot(Seurat_NK_Final, features = "percent.mt" )

Seurat_NK_Final = subset(Seurat_NK_Final, idents= c("0","1","2","3","4","7","8") )

DimPlot(Seurat_NK_Final, label=TRUE, label.size = 6)


#Have a look at their signatures
DefaultAssay(object = Seurat_NK_Final) <- "SCT"

All_Markers = FindAllMarkers(Seurat_NK_Final , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers at the transcriptional level


markers_overCl=All_Markers

All_Markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DotPlot(Seurat_NK_Final, features = unique(top10$gene) , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))



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


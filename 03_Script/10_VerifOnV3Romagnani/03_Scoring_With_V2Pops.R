

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/All_Markers_NK123.rds")

NUMBER_TOP_SCORING= 10

Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All

list_Top_Genes = list()
for (i in levels(top_All$cluster)){
  top_All %>%
    filter(cluster == i ) -> top_clust
  list_Top_Genes  =append(list_Top_Genes, list(top_clust$gene))
}

names(list_Top_Genes) = levels(top_All$cluster)
print(list_Top_Genes)

List_To_Use = lapply(list_Top_Genes , function(x) as.data.frame(x))

MONITORED_Markers = List_To_Use

PBMC= Merged_Seurat_Rescaled


for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}



#Plotting on a heatmap

PBMC$FirstClust = PBMC$seurat_clusters
PBMC$FirstClust  = droplevels(PBMC$FirstClust )
levels(PBMC$FirstClust) = c("Cluster0", "Cluster1" , "Cluster2" , "Cluster3" , "Cluster4" , "Cluster5", "Cluster6", "Cluster7")

table = PBMC@meta.data[,c( "FirstClust" , paste0(names(MONITORED_Markers),1) )]

df = table %>% 
  group_by(FirstClust) %>% 
  dplyr::summarise(across(everything(), list(mean)))

df= as.data.frame(df)
rownames(df)=df$FirstClust

df2= df[,-1]


pheatmap(as.matrix(df2), scale= "column", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_col = 15, fontsize_row = 15, show_rownames = TRUE, show_colnames = TRUE)

pheatmap(as.matrix(df2), scale= "column")

VlnPlot(PBMC , features = c("CXCR4", "XCL1", "XCL2", "SELL", "FCER1G"))
FeaturePlot(PBMC , features = c("CXCR4", "XCL1", "XCL2", "SELL", "FCER1G"))


DotPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), group.by = "seurat_clusters", cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(PBMC, features = "NK21", max.cutoff = 2.4)   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), split.by= "Majortype")   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC , features = paste0(names(MONITORED_Markers),"1"), group.by = "seurat_clusters", pt.size = 0) 

RidgePlot(PBMC , features = paste0(names(MONITORED_Markers),"1"), group.by = "celltype") 

FeaturePlot(PBMC, features= c("NKG7", "PRF1", "ITGB2", "SPON2", "GNLY" ))  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features= c("SELL", "CXCR4", "KLRB1", "FCGR3A", "SPON2", "FGFBP2", "FCER1G", "GZMH", "IL32" ))  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC, features= c("NKG7", "PRF1", "ITGB2", "SPON2", "GNLY" ), pt.size = 0, split.by = "Project")  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(PBMC, features= c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "CLIC3" ))  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(PBMC, features= c("KLRC2" ))  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))








j = 1
for (i in names(List_To_Use)){
  paste0("List_To_Use$" , i)
  #names(i) = names(List_To_Use)[j]
  j=j+1
}






colnames(List_To_Use$NK1A)
names(List_To_Use$NK1B)



DotPlot(Merged_Seurat_Rescaled, features = unique(top_All$gene) , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

DotPlot(Merged_Seurat_Rescaled, features = c("ACTB","SPON2",  "ACTG1", "PRF1", "PFN1") , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = c("CXCR4","KLRB1") , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = c("ZFP36","UTRN", "CLIC3", "CD160", "CD38") , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = c("ZFP36","UTRN", "CLIC3", "CD160", "CD38") , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

DotPlot(Merged_Seurat_Rescaled, features = MONITORED_Markers$NK1A , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = MONITORED_Markers$NK1B, cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = MONITORED_Markers$NK1C, cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = MONITORED_Markers$NK3, cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = MONITORED_Markers$NK2A, cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(Merged_Seurat_Rescaled, features = MONITORED_Markers$NK2B, cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

DotPlot(Merged_Seurat_Rescaled, features = unique(unlist(MONITORED_Markers)) , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))




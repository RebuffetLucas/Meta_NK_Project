#Scoring subcluster on the object with granularity 1.0 (=> 10 clusters in the end)

#For Dim

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP12_NK1.rds")


PBMC= Merged_Seurat_Rescaled

PBMC= subset(PBMC, idents= c("4", "5","1", "0") )

MONITORED_Markers = Markers_Seurat



for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#For Bright
PBMC= Merged_Seurat_Rescaled
PBMC= subset(PBMC, idents= c("7", "8") )

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_NK2.rds")


MONITORED_Markers = Markers_Seurat

for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#For NK3
PBMC= Merged_Seurat_Rescaled
PBMC= subset(PBMC, idents= c("0", "2", "3", "6", "9") )

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/All_Markers_NK3ABC.rds")


NUMBER_TOP_SCORING= 20

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


for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Marker based approach

FeaturePlot(PBMC, features= c("SELL","CXCR4", "KLRB1", "CD160","SPON2",  "CD3G", "CD3D", "GZMH", "IL32")) &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



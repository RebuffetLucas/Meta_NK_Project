#Scoring with our signatures

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/230626CVa_markers_vf1/resFirstClust_Seurat.RDS")

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/230626CVa_markers_vf1/resFirstClust_Seurat.RDS")
Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds")

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_NK2.rds")
Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP12_NK1.rds")

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_NK2.rds")

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/All_Markers_NK1ABC3.rds")




SUBSET = c("CD56brightCD16hi", "CD56dimCD16hi-c1-IL32"  , "CD56dimCD16hi-c2-CX3CR1", "CD56dimCD16hi-c3-ZNF90",  "CD56dimCD16hi-c4-NFKBIA" ,"CD56dimCD16hi-c8-KLRC2" ,
           "CD56brightCD16lo-c1-GZMH" ,"CD56brightCD16lo-c2-IL7R-RGS1lo", "CD56brightCD16lo-c3-CCL3", "CD56brightCD16lo-c4-IL7R")

PBMC = data_all

#Subseting
PBMC  =subset(PBMC, subset= celltype %in% SUBSET)

PBMC  =subset(PBMC, subset= Majortype == "CD56lowCD16high")
PBMC$celltype = droplevels(PBMC$celltype)

table(PBMC$celltype)


NUMBER_TOP_SCORING= 20

#Extracting top markers
Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All


#Building list of best markers
list_Top_Genes = list()
for (i in levels(top_All$cluster)){
  top_All %>%
    filter(cluster == i ) -> top_clust
  list_Top_Genes  =append(list_Top_Genes, list(top_clust$gene))
}

names(list_Top_Genes) = levels(top_All$cluster)
print(list_Top_Genes)


#Format that can be used ofr AddModuleScore
List_To_Use = lapply(list_Top_Genes , function(x) as.data.frame(x))
MONITORED_Markers = List_To_Use

for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}


#Several Plots

p1 = DimPlot(PBMC, group.by= "celltype" )
p2 = FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p1 + p2
DimPlot(PBMC ,group.by= "celltype", split.by= "Majortype" )
DimPlot(PBMC ,group.by= "celltype")


DotPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), group.by = "seurat_clusters", cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), group.by = "seurat_clusters", pt.size= 0)    & theme(axis.text.x = element_text(angle = 90))
FeaturePlot(PBMC, features = c("XCL1","XCL2", "LTB", "SELL"), split.by= "Majortype")   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(PBMC, features = c("KIR2DL1","KIR2DL3", "KIR3DL1", "KIR2DL5A"), split.by= "Majortype")   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC, features = c("KIR2DL1","KIR2DL3", "KIR3DL1", "KIR2DL5A"), group.by= "Majortype")   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = c("KIR2DL1","KIR2DL3", "KIR3DL1", "KIR2DL5A"), split.by= "Majortype")   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), split.by= "Majortype")   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


VlnPlot(PBMC , features = paste0(names(MONITORED_Markers),"1"), group.by = "celltype", pt.size = 0) 

RidgePlot(PBMC , features = paste0(names(MONITORED_Markers),"1"), group.by = "celltype") 



#Plotting on a heatmap


table = PBMC@meta.data[,c( "celltype" , paste0(names(MONITORED_Markers),1) )]

table = PBMC@meta.data[,c( "celltype" , paste0(names(MONITORED_Markers),1) )]


df = table %>% 
  group_by(celltype) %>% 
  dplyr::summarise(across(everything(), list(mean)))


df= as.data.frame(df)
rownames(df)=df$celltype
colnames(df)= gsub('.{3}$', '', colnames(df))
df2= df[,-1]


pheatmap(as.matrix(df2), scale= "column", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_col = 15, fontsize_row = 15, show_rownames = TRUE, show_colnames = TRUE)

pheatmap(as.matrix(df2), scale= "column")

VlnPlot(PBMC , features = c("CXCR4", "XCL1", "XCL2", "SELL", "FCER1G"), pt.size = 0, group.by = "celltype")



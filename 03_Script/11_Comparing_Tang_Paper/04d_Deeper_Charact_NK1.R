#Subset and find Variable features for NK 1


data_All_NK1 = subset(data_All, idents= c("4", "6", "0", "2", "1", "10"))

All_Markers = FindAllMarkers(data_All_NK1 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/Integrated_Seurat_Object/DEG_12clusters.rds")



#Look at markers like Constance
markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(data_All_NK1)))) {
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

topGenes2


topGenesb=unique(unlist(topGenes))

DotPlot(data_All_NK1, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

VlnPlot(data_All_NK1, features = c("nCount_RNA" , "nFeature_RNA" ), pt.size = 0) + theme(axis.text.x = element_text(angle = 90))


#Scoring_Intra_NK1

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/All_Markers_NK1ABC3.rds")


PBMC= data_All_NK1
#Malmberg= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg.rds")

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

DimPlot(PBMC)
DimPlot(PBMC ,group.by= "celltype", split.by= "Majortype" )
DimPlot(PBMC ,group.by=  "Project" )


DotPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), group.by = "seurat_clusters", cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1") , pt.size = 0)

DotPlot(PBMC, features=  c("ITGB2", "FCGR3A", "CD2", "NCR3", "FCER1G"), cols = "Spectral")

DotPlot(PBMC, features=  c("CXCR4", "JUNB", "CCL3", "CD38", "CLIC3", "CD160", "UTRN", "SYNE1", "ITGB2"), cols = "Spectral")


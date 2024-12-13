
#Import
PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")
PBMC_Meta_NK_V2= SetIdent(PBMC_Meta_NK_V2, value= "FirstClust")

print(DimPlot(PBMC_Meta_NK_V2, cols = palette))

#Subset NK1
NK1pop = subset(PBMC_Meta_NK_V2 , idents=c("NK1A", "NK1B", "NK1C"))

NK1andIntpop = subset(PBMC_Meta_NK_V2 , idents=c("NK1A", "NK1B", "NK1C", "NKint"))

#Diff exp Genes

All_Markers = FindAllMarkers(NK1pop , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


#saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/Integrated_Seurat_Object/DEG_12clusters.rds")
#Look at markers like Constance

All_Markers %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All



markers_overCl=top_All
topGenes=c()
for (cl in sort(unique(Idents(NK1pop)))) {
  print("Top genes for ")
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


topGenesb=unique(unlist(topGenes))


p10 = DotPlot(NK1pop, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))




All_Markers = FindAllMarkers(NK1andIntpop , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All



#saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/Integrated_Seurat_Object/DEG_12clusters.rds")
#Look at markers like Constance
markers_overCl=top_All
topGenes=c()
for (cl in sort(unique(Idents(NK1andIntpop)))) {
  print("Top genes for ")
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


topGenesb=unique(unlist(topGenes))

p11 = DotPlot(NK1andIntpop, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))


png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/21_FACS_Study/FACS_Study_Panels/Additional_BioInfo_Analysis_For_Markers/DotPlot_NK1ABC.png", width = 55, height = 15,  units = "cm", res=600 )
p10
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/21_FACS_Study/FACS_Study_Panels/Additional_BioInfo_Analysis_For_Markers/DotPlot_NK1ABCInt.png", width = 55, height = 15,  units = "cm", res=600 )
p11
dev.off()


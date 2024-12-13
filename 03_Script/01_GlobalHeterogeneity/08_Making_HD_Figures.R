#Making HD Figures for the paper

#def function

# +- std
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}



#Load objects
PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

PBMC  =SetIdent(PBMC, value = "FirstClust")



 
p1 = DimPlot(PBMC, label=TRUE, order = rev(ORDER_CLUST_LEGEND), cols = palette) + NoAxes()



png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/UMAP_6pops.png", width = 15, height = 10,  units = "cm", res=600 )
p1 
dev.off()

p_number = DimPlot(PBMC, label=TRUE, order = rev(ORDER_CLUST_LEGEND2), group.by = "SecondClust",cols = palette2) + NoAxes()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/UMAP_8pops_V2.png", width = 15, height = 10,  units = "cm", res=600 )
p_number
dev.off()

levels(PBMC$SecondClust) = c("2", "4", "0", "5", "6", "1", "8", "3")

palette_number = palette2
names(palette_number) = c("0", "5", "2", "4", "6", "1", "3", "8")

p_number2 = DimPlot(PBMC, label=TRUE, order = rev(ORDER_CLUST_LEGEND2), group.by = "SecondClust",cols = palette_number, ) + NoAxes()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/UMAP_8pops_V2_number.png", width = 15, height = 10,  units = "cm", res=600 )
p_number2
dev.off()


#Fig2a
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig2a.pdf",
  plot = p1,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 14,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



#Bar Diagram of Datasets, Chemistry , Orig.ident:

PBMC$seurat_clusters = PBMC$FirstClust
p5 = ggplot(PBMC@meta.data, aes(x=orig.ident, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette)
p7 = ggplot(PBMC@meta.data, aes(x=Dataset, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette)

p5 +  p7

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/BarGraph_Patients_6pops.png", width = 18, height = 18,  units = "cm", res=600 )
p5 
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/BarGraph_datasets_6pops.png", width = 12, height = 12,  units = "cm", res=600 )
p7
dev.off()

 
#Fig2b
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig2b.pdf",
  plot = p5,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 14,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



#Point plots for the proportion of 6 clusters in each patient
  #BoxPlot version
Clusters_Proportions = prop.table( table( PBMC$orig.ident , PBMC$seurat_clusters) , margin = 1)

Clusters_Proportions2 = t(Clusters_Proportions)
write.csv(Clusters_Proportions2, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/table_proportions_patients.csv", row.names=TRUE)


Clusters_Proportions = data.frame(Clusters_Proportions)

colnames(Clusters_Proportions) = c("orig.ident", "seurat_clusters", "Freq")

Clusters_Proportions$seurat_clusters = factor(Clusters_Proportions$seurat_clusters, levels= ORDER_CLUST_LEGEND)

p12 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_boxplot(outlier.shape =  NA, color = palette[ORDER_CLUST_LEGEND], lwd= 0.9, alpha= 0.1) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 3) +
  theme_classic() +  theme(text = element_text(size = 15))       


png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/MoustachePlot_Proportions_Patients_6pops.png", width = 20, height = 20,  units = "cm", res=600 )
p12
dev.off()


#Version BarPlot

p14 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", fill= "white", color= palette[ORDER_CLUST_LEGEND], size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.25, size=0.6) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() +  theme(text = element_text(size = 15))       

p14




png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/BarPlot_Proportions_Patients_6pops.png", width = 20, height = 20,  units = "cm", res=600 )
p14
dev.off()

#Do the same with NK1 , 2 and 3

PBMC$seurat_clusters2 = PBMC$seurat_clusters

levels(PBMC$seurat_clusters2) = c("NK1","NK1", "NK1", "NKint" ,"NK2", "NK3")

p11 = DimPlot(PBMC, group.by = "seurat_clusters2", cols= palette4, label= TRUE)

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/UMAP_4pops.png", width = 15, height = 10,  units = "cm", res=600 )
p11
dev.off()

PBMC$seurat_clusters2 = PBMC$seurat_clusters

levels(PBMC$seurat_clusters2) = c("NK1","NK1", "NK1", "NK1" ,"NK2", "NK3")

p10 = DimPlot(PBMC, group.by = "seurat_clusters2", cols= palette3, label= TRUE)

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/UMAP_3pops.png", width = 15, height = 10,  units = "cm", res=600 )
p10
dev.off()






Clusters_Proportions = prop.table( table( PBMC$orig.ident , PBMC$seurat_clusters2) , margin = 1)

Clusters_Proportions = data.frame(Clusters_Proportions)

colnames(Clusters_Proportions) = c("orig.ident", "seurat_clusters", "Freq")

Clusters_Proportions$seurat_clusters = factor(Clusters_Proportions$seurat_clusters, levels= ORDER_CLUST_LEGEND3)

p12 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_boxplot(outlier.shape =  NA, color = palette3[ORDER_CLUST_LEGEND3], lwd= 0.9, alpha= 0.1) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 3) +
  theme_classic() +  theme(text = element_text(size = 15))       

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/BarPlot_Proportions_Patients_3pops.png", width = 20, height = 20,  units = "cm", res=600 )
p12
dev.off()



#Point plots for the proportion of 8 clusters in each patient

Clusters_Proportions = prop.table( table( PBMC$orig.ident , PBMC$SecondClust) , margin = 1)

Clusters_Proportions = data.frame(Clusters_Proportions)

colnames(Clusters_Proportions) = c("orig.ident", "seurat_clusters", "Freq")

Clusters_Proportions$seurat_clusters = factor(Clusters_Proportions$seurat_clusters, levels= ORDER_CLUST_LEGEND2)

p20 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_boxplot(outlier.shape =  NA, color = palette2[ORDER_CLUST_LEGEND2], lwd= 0.9, alpha= 0.1) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 3) +
  theme_classic() +  theme(text = element_text(size = 15))       


png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/MoustachePlot_Proportions_Patients_8pops.png", width = 20, height = 20,  units = "cm", res=600 )
p20
dev.off()



#Unsupervised Clustering at low granularity
PBMC <- FindNeighbors(PBMC,  reduction = "harmony", dims= 1:30 )
PBMC <- FindClusters(PBMC, resolution = 0.1, random.seed = SEED )


DimPlot(PBMC)


#Scoring with CITEseq
PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

Markers_Seurat= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds") #Not ready for scoring

#write.csv(Markers_Seurat, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/20_Sup_Tables_csv_xlsx/DEG_CITEseq_NK1_NK2_NK3.csv", row.names=TRUE)
#write_xlsx(Markers_Seurat,"/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/20_Sup_Tables_csv_xlsx/DEG_CITEseq_NK1_NK2_NK3.xlsx")


#Extracting top markers
Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but not for plotind DotPlots
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


#Format that can be used by AddModuleScore
List_To_Use = lapply(list_Top_Genes , function(x) as.data.frame(x))
MONITORED_Markers = List_To_Use

for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

#VlnPlot and save
p3 = FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), max.cutoff = 1.5)   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
p3 = VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), pt.size = 0,  cols = palette) & stat_summary(fun.data=data_summary,color="black")

p4 = VlnPlot(PBMC, features = paste0(names(MONITORED_Markers),"1"), pt.size = 0,  cols = palette)

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/VlnPlot_Scoring_WNN_on35KUMAP.png", width = 30, height = 15,  units = "cm", res=600 )
p3 
dev.off()


#Fig2c
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig2c.pdf",
  plot = p3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#VlnPlot scoring NK3 in CMVneg vs CMVpos
#Extract only Romagnani Dataset
PBMC_Dataset4 = subset(PBMC, subset  = Dataset== "Dataset4")
PBMC_Dataset4 = subset(PBMC_Dataset4, subset  = FirstClust== "NK3")


#Separate CMVpos an neg
PBMC_Dataset4$orig.ident = droplevels(PBMC_Dataset4$orig.ident)
table(PBMC_Dataset4$orig.ident)

#Add CMV status
PBMC_Dataset4$CMVstatus = PBMC_Dataset4$orig.ident
levels(PBMC_Dataset4$CMVstatus) = c("CMVneg", "CMVneg", "CMVpos", "CMVpos", "CMVpos")
table(PBMC_Dataset4$CMVstatus)

p3 = VlnPlot(PBMC_Dataset4, features = "NK_31", group.by = "CMVstatus",pt.size = 0) & stat_summary(fun.data=data_summary,color="black")


p3

#Fig2d
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig2d.pdf",
  plot = p3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 15,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#Markers DotPlot
#Find Markers
PBMC@active.ident = factor(PBMC@active.ident, levels= ORDER_CLUST_LEGEND) 

PBMC= SetIdent(PBMC, value = "FirstClust")
All_Markers = FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


write.csv(All_Markers, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/20_Sup_Tables_csv_xlsx/DEG_scRNAseq_NK1ABC_NK2_NK3_NKint.csv", row.names=TRUE)
write_xlsx(All_Markers,"/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/20_Sup_Tables_csv_xlsx/DEG_scRNAseq_NK1ABC_NK2_NK3_NKint.xlsx")





#Look at markers like Constance
All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10


DotPlot(PBMC, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))
p4 = DotPlot(PBMC, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))

#Save that in HD
png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/DotPlot_Signatures_removeRPS_RPL_MT_RdBu.png", width = 55, height = 15,  units = "cm", res=600 )
p4
dev.off()

#Fig2e
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig2e.pdf",
  plot = p4,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 60,
  height = 12,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)





#For sup figure:

#DotPlot 8 clusters
test= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1.rds")
DimPlot(test)
levels(test@active.ident) = c("NK1C", "NK3A",  "NK1A", "NK3C" , "NK1B", "NKint", "NK2", "NK3B")
test@active.ident = factor(test@active.ident, levels= ORDER_CLUST_LEGEND2) 

p15 = DimPlot(test, cols = palette2, label = TRUE)+ NoAxes()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/UMAP_8clusters.png", width = 15, height = 10,  units = "cm", res=600 )
p15
dev.off()


All_Markers = FindAllMarkers(test , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE, cols= palette2)


All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10


DotPlot(test, features = unique(top10$gene), cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))
p5 =DotPlot(test, features = unique(top10$gene), cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/DotPlot_8clusters_Signatures_removeRPS_RPL_MT_RdBu.png", width = 55, height = 15,  units = "cm", res=600 )
p5
dev.off()

#DotPlot for the 3 pops of NK3
test$SecondClust = test@active.ident
test = subset(test, idents = c("NK3A", "NK3B", "NK3C") )

All_Markers = FindAllMarkers(test , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10


DotPlot(test, features = unique(top10$gene), cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
p5 = DotPlot(test, features =  unique(top10$gene) , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/DotPlot_NK3_3subpops_removeRPS_RPL_MT.png", width = 55, height = 15,  units = "cm", res=600 )
p5
dev.off()


#FeaturePlots


p16 = FeaturePlot(PBMC, features= c("GZMK", "XCL1",  "KLRC1", "SELL")) & NoAxes() &   theme( plot.title = element_text( face = "italic") )
p17 =FeaturePlot(PBMC, features= c("GZMB", "FCGR3A",  "SPON2", "FCER1G")) & NoAxes() &   theme( plot.title = element_text( face = "italic") )
p18 = FeaturePlot(PBMC, features= c("GZMH", "CCL5",  "KLRC2", "IL32")) & NoAxes() &   theme( plot.title = element_text( face = "italic") )

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/FeaturePlot_Markers_NK1.png", width = 20, height = 15,  units = "cm", res=600 )
p17
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/FeaturePlot_Markers_NK2.png", width = 20, height = 15,  units = "cm", res=600 )
p16
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/FeaturePlot_Markers_NK3.png", width = 20, height = 15,  units = "cm", res=600 )
p18
dev.off()



#Proportion of subset for Dataset4
PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

  #Extract only Romagnani Dataset
PBMC_Dataset4 = subset(PBMC, subset  = Dataset== "Dataset4")


#Separate CMVpos an neg
PBMC_Dataset4$orig.ident = droplevels(PBMC_Dataset4$orig.ident)
table(PBMC_Dataset4$orig.ident)

#Add CMV status
PBMC_Dataset4$CMVstatus = PBMC_Dataset4$orig.ident
levels(PBMC_Dataset4$CMVstatus) = c("CMVneg", "CMVneg", "CMVpos", "CMVpos", "CMVpos")
table(PBMC_Dataset4$CMVstatus)


#Point plots for the proportion of 8 clusters in each patient

Clusters_Proportions = prop.table( table( PBMC_Dataset4$orig.ident , PBMC_Dataset4$FirstClust) , margin = 1)
Clusters_Proportions = data.frame(Clusters_Proportions)

colnames(Clusters_Proportions) = c("orig.ident", "seurat_clusters", "Freq")

Clusters_Proportions$CMVstatus = Clusters_Proportions$orig.ident

levels(Clusters_Proportions$CMVstatus) = c("CMVneg", "CMVneg", "CMVpos", "CMVpos" ,"CMVpos")

table(Clusters_Proportions$CMVstatus , Clusters_Proportions$orig.ident)

Clusters_Proportions$seurat_clusters = factor(Clusters_Proportions$seurat_clusters, levels= ORDER_CLUST_LEGEND)

p14 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", fill= "white", color= palette[ORDER_CLUST_LEGEND], size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.25, size=0.6) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() +  theme(text = element_text(size = 15)) + facet_wrap(~ CMVstatus)

#This works better
p14 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", fill= "white", color= palette[ORDER_CLUST_LEGEND], size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.25, size=0.6) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() +  theme(text = element_text(size = 15))

p14




p14 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", fill= "white", color= palette[ORDER_CLUST_LEGEND], size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.25, size=0.6) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() +  theme(text = element_text(size = 15))

p14


p14 + facet_grid(~ CMVstatus)



p15 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, fill = CMVstatus, colour = orig.ident)) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), width = 0.25, size=0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  facet_grid(~ CMVstatus)

p15

p16 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, fill = CMVstatus, color = orig.ident)) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), width = 0.25, size=0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  facet_wrap(~ CMVstatus)

p16




p15 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, fill = CMVstatus)) + # Ajout de 'fill' pour la distinction CMVpos et CMVneg
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), color= palette[ORDER_CLUST_LEGEND], size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), width = 0.25, size=0.6) +
  geom_jitter(aes(color = orig.ident), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() + 
  theme(text = element_text(size = 15)) 



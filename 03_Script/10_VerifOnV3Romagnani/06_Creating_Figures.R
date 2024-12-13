
#Create figures

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/3pops_Object.rds")

NUMBER_TOP_SCORING = 20

p6 = FeaturePlot(PBMC, features= c("GZMK","SELL", "XCL1",
                                   "GZMB" , "FCGR3A", "SPON2", 
                                   "GZMH", "CD52",  "IL32" ))  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

p7 = FeaturePlot(PBMC, features= c("GZMK","SELL", "XCL1",
                                   "GZMB" , "FCGR3A", "SPON2", 
                                   "GZMH", "CD52",  "IL32" ))  & NoAxes() &   theme( plot.title = element_text( face = "italic") )




png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/FeaturePlot_NKMarkers.png", width = 25, height = 20,  units = "cm", res=600 )
p6
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/FeaturePlot_NKMarkers_Violet.png", width = 25, height = 20,  units = "cm", res=600 )
p7
dev.off()



#Markers

All_Markers = FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All


All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10


DotPlot(PBMC, features = unique(top10$gene), cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))
p4 = DotPlot(PBMC, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/DotPlot_V3_Signature_RdBu.png", width = 35, height = 10,  units = "cm", res=600 )
p4
dev.off()

Genes_Plot= unique(top10$gene)[c(2,5,17,18,20,    21,23,25,28,30,   41,42,43,44,45 )]

Genes_Plot= unique(top10$gene)


CITE_seq_Signature = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds")

CITE_seq_Signature %>%
  #filter(pct.1>0.5) %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_CITE_seq


top_CITE_seq  %>%
  filter(cluster == "NK_1") -> top_CITE_seq_NK1

top_CITE_seq  %>%
  filter(cluster == "NK_2") -> top_CITE_seq_NK2

top_CITE_seq  %>%
  filter(cluster == "NK_3") -> top_CITE_seq_NK3

#genes_to_plot = Reduce(union, list(top_CITE_seq_NK1$gene, top_CITE_seq_NK2$gene, top_CITE_seq_NK3$gene))


genes_to_plot = Reduce(union, list(top_CITE_seq_NK1$gene, top_CITE_seq_NK3$gene , top_CITE_seq_NK2$gene))
 
genes_to_plot2 = Reduce(intersect, list(genes_to_plot, top10$gene))



p8 = DotPlot(PBMC, features =  genes_to_plot2 , cols = "RdBu") + coord_flip()
p9 = DotPlot(PBMC, features =  genes_to_plot2 , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) 

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/DotPlot_V3_Shared_Signature_Vertical_RdBu.png", width = 12, height = 35,  units = "cm", res=600 )
p8
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/DotPlot_V3_Shared_Signature_Horizontal_RdBu.png", width = 35, height = 8,  units = "cm", res=600 )
p9
dev.off()



# Scoring with our signatures

#Scoring with 35K UMAP
Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds") #Already ready for scoring

MONITORED_Markers = Markers_Seurat

for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

p1 = DimPlot(PBMC )  & NoAxes()
p2 = FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  & NoAxes()


p2 + p1

#Save Figures

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/Scoring_V3_With_V2_3pops.png", width = 15, height = 12,  units = "cm", res=600 )
p2 
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/UMAP_V3_3pops.png", width = 15, height = 12,  units = "cm", res=600 )
p1 
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/DotPlot_V3_shared_Markers.png", width = 15, height = 12,  units = "cm", res=600 )
p8 
dev.off()





#Scoring with CITEseq
PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/3pops_Object.rds")

Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds") #Already ready for scoring

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

#Plot and save
p3 = FeaturePlot(PBMC, features = paste0(names(MONITORED_Markers),"1"))   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  & NoAxes()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/14_Verif_V3_Romagnani/Final_Figures/UMAP_V3_scoredWithWNN_3pops.png", width = 15, height = 12,  units = "cm", res=600 )
p3 
dev.off()


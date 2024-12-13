#Scripts to make he figures that were finally retained

#Parameters
NUMBER_TOP_SCORING= 20

#Loading the objects

PBMC_Blood = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5seurat")
PBMC_All = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5seurat")

#Preparing the objects
SUBSET = c("CD56brightCD16hi", "CD56dimCD16hi-c1-IL32"  , "CD56dimCD16hi-c2-CX3CR1", "CD56dimCD16hi-c3-ZNF90",  "CD56dimCD16hi-c4-NFKBIA" ,"CD56dimCD16hi-c8-KLRC2" ,
           "CD56brightCD16lo-c1-GZMH" ,"CD56brightCD16lo-c2-IL7R-RGS1lo", "CD56brightCD16lo-c3-CCL3", "CD56brightCD16lo-c4-IL7R")
table(PBMC_Blood$celltype )
levels(PBMC_Blood$celltype) = levels(PBMC_All$celltype)

PBMC = PBMC_Blood


table(PBMC$celltype)


#Subseting to keep only the >1% clusters

PBMC  =subset(PBMC, subset= celltype %in% SUBSET)
PBMC$celltype = droplevels(PBMC$celltype)

table(PBMC$celltype)


#Scoring with 35K UMAP
Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds") #Already ready for scoring

MONITORED_Markers = Markers_Seurat

for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

p1 = DimPlot(PBMC, group.by= "celltype" ) + NoAxes()
p2 = FeaturePlot(PBMC, features = "NK11") + NoAxes() + ggtitle("NK1")  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 = FeaturePlot(PBMC, features = "NK21") + NoAxes()  + ggtitle("NK2") &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p4 = FeaturePlot(PBMC, features = "NK31") + NoAxes()  + ggtitle("NK3") &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p1 + p2 + p3 + p4

p = p1 + p2 + p3 + p4

#Save Figures

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/FinalFigures/Scoring_With_V2_3pops/Scoring_Tang_With_V2_3pops.png", width = 25, height = 15,  units = "cm", res=600 )
p
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/FinalFigures/Tang_Blood_UMAP.png", width = 22, height = 13,  units = "cm", res=600 )
p1 
dev.off()


p5 = RidgePlot(PBMC , features = paste0(names(MONITORED_Markers),"1"), group.by = "celltype") 

p6 = RidgePlot(PBMC, features = "NK11", group.by = "celltype") + ggtitle("NK1")  + theme(axis.text = element_text(size = 12))     
p7 = RidgePlot(PBMC, features = "NK21", group.by = "celltype")   + ggtitle("NK2")  + theme(axis.text = element_text(size = 12)) 
p8 = RidgePlot(PBMC, features = "NK31", group.by = "celltype")  + ggtitle("NK3") + theme(axis.text = element_text(size = 12))


p6 + p7 + p8
p5 = p6 + p7 +p8

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/FinalFigures/Scoring_With_V2_3pops/Scoring_Tang_With_V2_3pops_RidgePlot.png", width = 25, height = 30,  units = "cm", res=600 )
p5 
dev.off()


#Scoring with CITEseq
PBMC = PBMC_Blood
table(PBMC$celltype)


#Subseting to keep only the >1% clusters

PBMC  =subset(PBMC, subset= celltype %in% SUBSET)
PBMC$celltype = droplevels(PBMC$celltype)
table(PBMC$celltype)


#Scoring with 35K UMAP
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
p1 = DimPlot(PBMC, group.by= "celltype" ) + NoAxes()
p2 = FeaturePlot(PBMC, features = "NK_11") + NoAxes() + ggtitle("NK1")  &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 = FeaturePlot(PBMC, features = "NK_21") + NoAxes()  + ggtitle("NK2") &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p4 = FeaturePlot(PBMC, features = "NK_31") + NoAxes()  + ggtitle("NK3") &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p1 + p2 + p3 + p4

p = p1 + p2 + p3 + p4

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/FinalFigures/Scoring_With_WNN_3pops/Scoring_Tang_With_WNN_3pops.png", width = 25, height = 15,  units = "cm", res=600 )
p 
dev.off()

p6 = RidgePlot(PBMC, features = "NK_11", group.by = "celltype") + ggtitle("NK1")  + theme(axis.text = element_text(size = 12))     
p7 = RidgePlot(PBMC, features = "NK_21", group.by = "celltype")   + ggtitle("NK2")  + theme(axis.text = element_text(size = 12)) 
p8 = RidgePlot(PBMC, features = "NK_31", group.by = "celltype")  + ggtitle("NK3") + theme(axis.text = element_text(size = 12))


p6 + p7 + p8
p5 = p6 + p7 +p8

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/FinalFigures/Scoring_With_WNN_3pops/Scoring_Tang_With_WNN_3pops_RidgePlot.png", width = 25, height = 30,  units = "cm", res=600 )
p5
dev.off()





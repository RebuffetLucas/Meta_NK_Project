
PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

Markers_3pops_CITEseq = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds")

NUMBER_TOP_SCORING = 20


#Extract only Romagnani Dataset
PBMC_Dataset4 = subset(PBMC, subset  = Dataset== "Dataset4")

#Separate CMVpos an neg

PBMC_Dataset4$orig.ident = droplevels(PBMC_Dataset4$orig.ident)
table(PBMC_Dataset4$orig.ident)

PBMC_Dataset4$CMVstatus = PBMC_Dataset4$orig.ident
levels(PBMC_Dataset4$CMVstatus) = c("CMVneg", "CMVneg", "CMVpos", "CMVpos", "CMVpos")

table(PBMC_Dataset4$CMVstatus)

#Extract best genes from CITEseq
Markers_Seurat= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds") #Not Already ready for scoring

#Extracting top markers
Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
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

names(list_Top_Genes)= c("NK_1", "NK_3", "NK_2")

  #Format that can be used by AddModuleScore
List_To_Use = lapply(list_Top_Genes , function(x) as.data.frame(x))
MONITORED_Markers = List_To_Use

for (i in names(MONITORED_Markers)){
  PBMC_Dataset4 = AddModuleScore(PBMC_Dataset4, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

#Score and Plot

  #Dataset4 all cells
p1 = VlnPlot(PBMC_Dataset4, features= c("NK_11", "NK_21","NK_31"), group.by = "CMVstatus", pt.size = 0)

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/VlnPlot_CMVposvsNeg_Dataset4_Allsubsets_ScoredwithCITEseqSignatures.png", width = 18, height = 12,  units = "cm", res=600 )
p1
dev.off()

  #Dataset4 Only NK3 cells

PBMC_Dataset4_NK3 = subset(PBMC_Dataset4, subset  =  FirstClust== "NK3")


p2 = VlnPlot(PBMC_Dataset4_NK3, features= c("NK_11", "NK_21","NK_31"), group.by = "CMVstatus", pt.size = 0)

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/VlnPlot_CMVposvsNeg_Dataset4_NK3subset_ScoredwithCITEseqSignatures.png", width = 18, height = 12,  units = "cm", res=600 )
p2
dev.off()




PBMC_Dataset4_NK1 = subset(PBMC_Dataset4, subset  =  FirstClust== c("NK1A", "NK1B", "NK1C"))

table(PBMC_Dataset4_NK1$FirstClust)

p3 = VlnPlot(PBMC_Dataset4_NK1, features= c("NK_11", "NK_21","NK_31"), group.by = "CMVstatus", pt.size = 0)

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/VlnPlot_CMVposvsNeg_Dataset4_NK1subset_ScoredwithCITEseqSignatures.png", width = 18, height = 12,  units = "cm", res=600 )
p3
dev.off()



p4 = VlnPlot(PBMC_Dataset4_NK1, features= c("NK_11", "NK_21","NK_31"), group.by = "FirstClust", pt.size = 0, split.by = "CMVstatus")

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/VlnPlot_CMVposvsNeg_Dataset4_NK1subset_Split1A1B1C_ScoredwithCITEseqSignatures.png", width = 18, height = 12,  units = "cm", res=600 )
p4
dev.off()







#Exploring alternatives


### Analysis of SCENIC Output ####

#Load packages
#For analysis
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
#For some of the plots
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(loomR)
library(dplyr)



number_regulons_showtop = 30
Z_TRESHOLD = 0.0000001
RegulonsDvt= c("EOMES(+)","STAT5A(+)","STAT5B(+)","TBX21(+)","PRDM1(+)","GATA3(+)","SMAD4(+)","FOXO1(+)", "NFIL3(+)", "ETS1(+)", "FOXO1(+)", "TCF7(+)", "ZEB1(+)", "MYC(+)", "IRF8", "RUNX3")

RegulonsOfInterest= c("EOMES(+)","STAT5A(+)","STAT5B(+)","TBX21(+)","PRDM1(+)","GATA3(+)","SMAD4(+)","FOXO1(+)", "NFIL3(+)", "ETS1(+)", "FOXO1(+)", "TCF7(+)", "ZEB1(+)", "MYC(+)", "IRF8(+)",
                      "RUNX3(+)", "RARA(+)","RUNX2(+)", "ZEB1(+)", "KLF6(+)", "ASCL2(+)", "NFIL3(+)", "IKZF1(+)", "TBX21(+)", "KLF2(+)","REL(+)","CREM(+)","FOXO4(+)","FOXO3(+)","KLF3(+)","BCL3(+)") 

RegulonsOfInterest2 = gsub('.{3}$', '', RegulonsOfInterest)



############Extract the names of True Human TF ##############################"""

library(readxl)
Can_TF= read_excel("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/List_Canonical_human_TF.xlsx")

# Assuming your tibble is named Can_TF
filtered_tibble <- Can_TF %>%
  select(ID, Name, DBD, TF, `TF assessment`) %>%
  filter(TF == "Yes" & `TF assessment` %in% c("Inferred motif", "Known motif", "Likely to be sequence specific TF"))

# View the filtered tibble
head(filtered_tibble)
table(filtered_tibble$`TF assessment`)

Name_Canonical_Human_TF = filtered_tibble$Name

######################################################################################




#For VF2:


MetaNKinfo=readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK_Info.rds")

MetaNKinfo$FirstClust
table(MetaNKinfo$FirstClust)

levels(MetaNKinfo$FirstClust) = c("NK1C", "NK3", "NK1A", "NK1B", "NKint", "NK2")


vsnDir <- "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/AllNK"

scenicLoomPath <- paste0(vsnDir, "/auc_mtx.loom")
motifEnrichmentFile <- paste0(vsnDir, "/expr_mat.adjacencies.tsv")
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)




#Loading the initial loom file
loom <- open_loom(scenicLoomPath, mode="r+")

regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

#Add the results of SCENIC analysis in the loom
#read info from loom file
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) #Better with log normalization
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAUCTresholds = get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
#cellClusters <- get_clusterings(loom)    #Indisponible sur cette analyse

close_loom(loom)

#Check before analysis

length(regulons)
head(names(regulons))
regulonAUC



#Load the motif enrichment results

motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=0)
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")



##Scoring the network activity

#Regulators for known cell types or clusters
#Average Regulon Activity per cluster


cellClusterPrecise=MetaNKinfo[,"FirstClust", drop=FALSE]


length(colnames(PBMC))
length(colnames(regulonAUC))
length(rownames(motifEnrichment))



#For now only:
Cells_to_keep  = intersect(rownames(MetaNKinfo),  colnames(regulonAUC) )
cellClusterPrecise=MetaNKinfo[Cells_to_keep,"FirstClust", drop=FALSE]
selectedResolution <- "FirstClust" # select resolution

regulonAUC = regulonAUC[,Cells_to_keep] 
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusterPrecise[Cells_to_keep, selectedResolution])


## ANALYSIS FOR PRECISE CLUSTERS

selectedResolution <- "FirstClust" # select resolution


# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusterPrecise), cellClusterPrecise[,selectedResolution]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

#Remove the "+"
rownames(regulonAUC) = gsub('.{3}$', '', rownames(regulonAUC))
######################## SUBSET ONLY REGULONS CORRESPONDING TO REAL HUMAN TF #####################################"

TF_To_Keep = intersect( rownames(regulonAUC) , Name_Canonical_Human_TF)
length(TF_To_Keep)


regulonAUC <- regulonAUC[TF_To_Keep,]

# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#Horizontal heatmap with all regulons
# Scale expression:
  #regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_byCellType_Scaled <- scale(t(regulonActivity_byCellType), center = T, scale=T)
#colnames(regulonActivity_byCellType_Scaled) = gsub('.{3}$', '', colnames(regulonActivity_byCellType_Scaled))



# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size



#Easy basic
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", column_names_gp=grid::gpar(fontsize=7), 
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later



png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/Final_Figures_V2/heatmap_AllREGULONS_Horizontal_TrueTF.png", width = 50, height = 10,  units = "cm", res=600 )
draw(hm)
dev.off()

#Same_By_Column:

regulonActivity_byCellType_Scaled_colones = t(regulonActivity_byCellType_Scaled)
hm2 <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled_colones, name="Regulon activity", column_names_gp=grid::gpar(fontsize=15), 
                                   row_names_gp=grid::gpar(fontsize=7))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/Final_Figures_V2/heatmap_AllREGULONS_Vertical_TrueTF.png", width = 10, height = 40,  units = "cm", res=600 )
draw(hm2)
dev.off()





#Horizontal heatmap with only regulons of interest genes ploted

test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""


hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=7), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/heatmap_AllREGULONS_Horizontal_OnlygenesOfInterest_TrueTF.png", width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()





#Horizontal heatmap including only the regulons of interest

test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""


hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=7), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/heatmap_AllREGULONS_Horizontal_OnlygenesOfInterest_TrueTF.png", width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()


#Plot for regulons of interest

RegulonsOfInterest2 = intersect(colnames(regulonActivity_byCellType_Scaled), RegulonsOfInterest2)


DvtRegulonsplot= regulonActivity_byCellType_Scaled[,unique(as.character(RegulonsOfInterest2))]

pheatmap::pheatmap(DvtRegulonsplot, fontsize_row=12, fontsize_col = 10, angle_col= 90,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey") 

hm2 <- draw(ComplexHeatmap::Heatmap(DvtRegulonsplot, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=8), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=10))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/Final_Figures_V2/heatmap_RegulonsOfInterest_Horizontal_TrueTF.png", width = 15, height = 10,  units = "cm", res=600 )
draw(hm2)
dev.off()


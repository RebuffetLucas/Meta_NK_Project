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


number_regulons_RSS_showtop= 10
number_regulons_showtop = 30
Z_TRESHOLD = 0.0000001
RegulonsDvt= c("EOMES(+)","STAT5A(+)","STAT5B(+)","TBX21(+)","PRDM1(+)","GATA3(+)","SMAD4(+)","FOXO1(+)", "NFIL3(+)", "ETS1(+)", "FOXO1(+)", "TCF7(+)", "ZEB1(+)", "MYC(+)", "IRF8", "RUNX3")

RegulonsOfInterest= c("EOMES(+)","STAT5A(+)","STAT5B(+)","TBX21(+)","PRDM1(+)","GATA3(+)","SMAD4(+)","FOXO1(+)", "NFIL3(+)", "ETS1(+)", "FOXO1(+)", "TCF7(+)", "ZEB1(+)", "MYC(+)", "IRF8(+)",
                      "RUNX3(+)", "RARA(+)","RUNX2(+)", "ZEB1(+)", "KLF6(+)", "ASCL2(+)", "NFIL3(+)", "IKZF1(+)", "TBX21(+)", "KLF2(+)","REL(+)","CREM(+)","FOXO4(+)","FOXO3(+)","KLF3(+)","BCL3(+)") 

RegulonsOfInterest2 = gsub('.{3}$', '', RegulonsOfInterest)

#For VF2:


MetaNKinfo=readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK_Info.rds")

MetaNKinfo$FirstClust
table(MetaNKinfo$FirstClust)

levels(MetaNKinfo$FirstClust) = c("NK1C", "NK3", "NK1A", "NK1B", "NKint", "NK2")

table(MetaNKinfo$FirstClust)


vsnDir <- "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/AllNK"

scenicLoomPath <- paste0(vsnDir, "/auc_mtx.loom")
motifEnrichmentFile <- paste0(vsnDir, "/expr_mat.adjacencies.tsv")
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)




#Loading the initial loom file
loom <- open_loom(scenicLoomPath, mode="r+")

regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

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


# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
#Horizontal heatmap with all regulons
# Scale expression:
  #regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_byCellType_Scaled <- scale(t(regulonActivity_byCellType), center = T, scale=T)
colnames(regulonActivity_byCellType_Scaled) = gsub('.{3}$', '', colnames(regulonActivity_byCellType_Scaled))




# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size



#Easy basic
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", column_names_gp=grid::gpar(fontsize=5), 
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later



png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/Final_Figures_V2/heatmap_AllREGULONS_Horizontal.png", width = 50, height = 10,  units = "cm", res=600 )
draw(hm)
dev.off()

#Same_By_Column:

regulonActivity_byCellType_Scaled_colones = t(regulonActivity_byCellType_Scaled)
hm2 <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled_colones, name="Regulon activity", column_names_gp=grid::gpar(fontsize=15), 
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/Final_Figures_V2/heatmap_AllREGULONS_Vertical.png", width = 10, height = 40,  units = "cm", res=600 )
draw(hm2)
dev.off()





#Horizontal heatmap with only regulons of interest genes ploted

test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""


hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=7), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/heatmap_AllREGULONS_Horizontal_OnlygenesOfInterest.png", width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()





#Horizontal heatmap including only the regulons of interest

test = regulonActivity_byCellType_Scaled

colnames(test)[!(colnames(test) %in% RegulonsOfInterest2)] = ""


hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=7), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/heatmap_AllREGULONS_Horizontal_OnlygenesOfInterest.png", width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()


#Plot for regulons of interest



RegulonsOfInterest2 = intersect(colnames(regulonActivity_byCellType_Scaled), RegulonsOfInterest2)


DvtRegulonsplot= regulonActivity_byCellType_Scaled[,unique(as.character(RegulonsOfInterest2))]

pheatmap::pheatmap(DvtRegulonsplot, fontsize_row=12, fontsize_col = 10, angle_col= 90,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey") 

hm2 <- draw(ComplexHeatmap::Heatmap(DvtRegulonsplot, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=5), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=10))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/Final_Figures_V2/heatmap_RegulonsOfInterest_Horizontal.png", width = 15, height = 10,  units = "cm", res=600 )
draw(hm2)
dev.off()





#Same with AnnoMArk:

test = regulonActivity_byCellType_Scaled
colnames(test) = gsub('.{3}$', '', colnames(test))

indices = which(colnames(test) %in% RegulonsOfInterest2)

labels_of_interest = colnames(test)[indices]

ha = columnAnnotation(foo = anno_mark(at = indices, 
                                   labels = labels_of_interest))

hm <- draw(ComplexHeatmap::Heatmap(test, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=0), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1), bottom_annotation = ha, column_names_side = "bottom", 
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size



#See the exact values

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

viewTable(topRegulators, options = list(pageLength = 10))

#Cell-type specific regulators

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusterPrecise[colnames(regulonAUC), selectedResolution])

## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss, zThreshold= Z_TRESHOLD)
plotly::ggplotly(rssPlot$plot)


options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = "NK1A", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "NK1B", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "NK1C", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "NK2", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "NK3", n=number_regulons_RSS_showtop) # cluster ID
plotRSS_oneSet(rss, setName = "NKint", n=number_regulons_RSS_showtop) # cluster ID


number_regulons_showtop = 10
#Extract the best Regulons and do a heatmap
rssPlot$df %>%
  dplyr::group_by(cellType) %>%
  top_n(n = number_regulons_showtop, wt = RSS) -> top10Regulons


regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]

unique(as.character(top10Regulons$Topic))

pheatmap::pheatmap(topregulonsplot, fontsize_row=10, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons


topregulonsplot2 = t(topregulonsplot)

colnames(topregulonsplot2) = gsub('.{3}$', '', colnames(topregulonsplot2))

hm <- draw(ComplexHeatmap::Heatmap(topregulonsplot2, name="Regulon activity",  column_names_gp=grid::gpar(fontsize=9), rect_gp = gpar(col = "white", lwd = 0.4), border_gp = gpar(col = "darkgrey", lty = 1),
                                   row_names_gp=grid::gpar(fontsize=15))) # row font size

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/Final_Figures/heatmap_AllREGULONS_Horizontal_OnlygenesOfInterest_TOP10.png", width = 50, height = 15,  units = "cm", res=600 )
draw(hm)
dev.off()





rssPlot$df %>%
  top_n(n = number_Best_regulons, wt = RSS) -> BestRegulons

unique(as.character(BestRegulons$Topic))

bestregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(BestRegulons$Topic)),]

pheatmap::pheatmap(bestregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons



#Look at all with RSS>0.01
rssPlot$df -> top10Regulons

topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]


pheatmap::pheatmap(topregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons


#Plot for dvt of NK cells


RegulonsDvt= c("EOMES(+)","STAT5A(+)","STAT5B(+)","TBX21(+)","PRDM1(+)","GATA3(+)","SMAD4(+)","FOXO1(+)", "NFIL3(+)", "ETS1(+)", "FOXO1(+)", "TCF7(+)", "ZEB1(+)", "MYC(+)", "IRF8(+)", "RUNX3(+)")

RegulonsDvt = intersect(rownames(regulonActivity_byCellType_Scaled), RegulonsDvt)


DvtRegulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(RegulonsDvt)),]

pheatmap::pheatmap(DvtRegulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons




#### Vizualisation with UMAP etc... => Has to done, necessit embedings ! 


#Cell states based on the GRN activity

# List of embeddings available:
cat(names(embeddings), sep="\n")

# Overview of the embeddings (see below for details)

#Plot using AUCell_plotTSNE
regulonsToPlot <- "STAT5B(+)"


#Plot cells present in embeddings and regulonAUC:
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "TCF7"
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))




# Plot cells present in both embeddings an regulonAUC
Embeddings_to_Plot = embeddings[["UMAP"]][colnames(regulonAUC),]

Embeddings_to_Plot[,"_X"] = -Embeddings_to_Plot[,"_X"]
Embeddings_to_Plot[,"_Y"] = -Embeddings_to_Plot[,"_Y"]

colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(Embeddings_to_Plot, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 


str(exprMat_log)


length(intersect(rownames(embeddings[["UMAP"]][colnames(regulonAUC),]) , colnames(regulonAUC) ))


  
  
  
#Plot cells present in embeddings and regulonAUC:
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "TCF7(+)"
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))




#Plot all cells in UMAP
colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "TCF7"
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))


#Plot of Bertrand # A retravailler car pas ouf du tout
df<-data.frame(UMAP_1=embeddings[["UMAP"]][,1],UMAP_2=embeddings[["UMAP"]][,2])
df$regulon<-regulonAUC@assays@data$AUC[regulonsToPlot,]
ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=regulon))+geom_point()+ scale_color_gradient2(low = "blue",high = "darkred",mid = "white",midpoint = max(df$regulon)/2)+theme_classic()

head(df)



##### The network in details: TFs, targets and motifs

length(regulons)
sum(lengths(regulons)>=10)


viewTable(cbind(nGenes=lengths(regulons)), options=list(pageLength=10))


# Check if a specific gene has a regulon


grep("EOMES", names(regulons), value=T) # paste0("^","EBF1","_")
grep("FOXP1", names(regulons), value=T) # paste0("^","EBF1","_")


#Look at the potential target genes

regulons[["EOMES(+)"]]
regulons[["TGIF2(+)"]]
regulons[["ZNF628(+)"]]
regulons[["RORA(+)"]]

regulons[["FOXP1(+)"]]
regulons[["KMT2A(+)"]]
regulons[["STAT4(+)"]]

regulons[["REL(+)"]]
regulons[["CREM(+)"]]


#Potential regulators for a given gene
gene <- "CD69"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "RGS1"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "IER2"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]

gene <- "IER5"
names(regulons)[which(sapply(regulons, function(x) gene %in% x))]


#dO IT WITH MULTIPLE GENES

dim(regulons_incidMat)

genes <- c("CD3E", "KLRF1", "FCGR3A") 
incidMat_subset <- regulons_incidMat[,genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset)>0,]

incidMat_subset


# Motifs supporting the regulons

tableSubset <- motifEnrichment[TF=="ZNF683"]
viewMotifs(tableSubset, colsToShow = c("logo", "NES", "TF" ,"Annotation"), options=list(pageLength=5))

head(tableSubset)

#Regulon targets and motifs

regulons[ c("EOMES(+)", "TBX21(+)")]


#Regulators for clusters of known cell types
#Average regulon activity by cluster

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)



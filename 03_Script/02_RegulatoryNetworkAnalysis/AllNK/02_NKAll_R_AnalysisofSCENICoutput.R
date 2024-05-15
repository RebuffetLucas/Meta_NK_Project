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



#For VF2:


MetaNKinfo=readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK_Info.rds")

MetaNKinfo$FirstClust

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

# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)


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
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

number_regulons = 20
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = "iNK1", n=number_regulons) # cluster ID
plotRSS_oneSet(rss, setName = "iNK2", n=number_regulons) # cluster ID
plotRSS_oneSet(rss, setName = "mNK1", n=number_regulons) # cluster ID
plotRSS_oneSet(rss, setName = "mNK2", n=number_regulons) # cluster ID
plotRSS_oneSet(rss, setName = "mNK3", n=number_regulons) # cluster ID



#Extract the best Regulons and do a heatmap
rssPlot$df %>%
  dplyr::group_by(cellType) %>%
  top_n(n = number_regulons, wt = RSS) -> top10Regulons



topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]

pheatmap::pheatmap(topregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons


#Look at all with RSS>0.01
rssPlot$df -> top10Regulons

topregulonsplot= regulonActivity_byCellType_Scaled[unique(as.character(top10Regulons$Topic)),]


pheatmap::pheatmap(topregulonsplot, fontsize_row=12, fontsize_col = 20,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color= "grey")  #Average expression regulon activity for top10Regulons


#Plot for dvt of NK cells


RegulonsDvt= c("EOMES(+)","STAT5A(+)","STAT5B(+)","TBX21(+)","PRDM1(+)","GATA3(+)","SMAD4(+)","FOXO1(+)", "NFIL3(+)", "ETS1(+)", "FOXO1(+)", "TCF7(+)", "ZEB1(+)", "MYC(+)", "IRF8", "RUNX3")

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
regulonsToPlot <- "TCF7(+)"

colorpalet=brewer.pal(n=8, name="RdBu")
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols =  colorpalet, offColor = "lightgray") 

regulonsToPlot <- "EOMES(+)"
AUCell::AUCell_plotTSNE(embeddings[["UMAP"]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, exprCols = c("goldenrod1", "darkorange", "brown"))


#Plot of Bertrand # A retravailler car pas ouf du tout
df<-data.frame(UMAP_1=embeddings[["UMAP"]][,1],UMAP_2=embeddings[["UMAP"]][,2])
df$regulon<-regulonAUC@assays@data$AUC[regulonsToPlot,]
ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=regulon))+geom_point()+scale_color_gradient2(low = "blue",high = "darkred",mid = "white",midpoint = max(df$regulon)/2)+theme_classic()





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


    
# NK vs NK2 vs NK3

#Perform GO analysis on a given seurat Object

#Part 1:
#  The aim of this script is to retrieve Functionnal Annotation Datasets
#  & Format both these datasets, DEG lists and Background Lists for ORA


library(clusterProfiler)
library(rlist)



PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")
UNIVERSE_SymbolID = rownames(PBMC@assays[["RNA"]])

THR_adj_pval = 0.05
THR_log2FC = 0

FX_TO_RUN = "groupGO" # One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .


Markers= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/DEG/DEG_NK1vsNK2vsNK3/MarkersNK1.2.3.rds")


# Search for ENTREZID equivalent ID
##### ADD ENTREZID TO ORA DF #####
# Rename Marker column to SYMBOL
ORA_DF <- Markers %>% dplyr::select(cluster,gene) %>% dplyr::rename(SYMBOL=gene)


# Conversion Gene Symbol to ENTREID for some functional enrichment analysis (such enrichKEGG)
genes_EntrezID<-AnnotationDbi::select(org.Hs.eg.db, ORA_DF$SYMBOL, 'ENTREZID','SYMBOL')
genes_EntrezID<-genes_EntrezID[which(genes_EntrezID$ENTREZID != is.na(genes_EntrezID$ENTREZID)),]
ORA_DF<-ORA_DF %>% filter(SYMBOL %in% genes_EntrezID$SYMBOL)
ORA_DF<-inner_join(ORA_DF,genes_EntrezID, by="SYMBOL") %>% dplyr::distinct()

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# CLUSTER geneList retrieval
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr geneList_retrieval

# Create an empty list of list for DLBCLs metaclusters marker genes 
ORA_DF_clusters <- vector(mode = "list", length = length(unique(ORA_DF$Metacluster)))
names(ORA_DF_clusters) <- unique(ORA_DF$Metacluster)

# Fill ORA_DF_clusters lists for each cutree with metacluster marker genes
for(j in seq_along(ORA_DF_clusters)){
  ORA_DF_clusters[[j]] <- ORA_DF %>% filter(Metacluster == paste0("MCl.",j))
}


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# pathways loading
# ---------------------------------------------------------------
# ---------------------------------------------------------------
## @knitr pathways_loading

### Read databases and create data.frames

#staudt_signatures <- as.data.frame(read_excel(STAUDT_REFERENCE_FILE))
#staudt_signatures <- staudt_signatures[,c(2,10)]
#colnames(staudt_signatures) <- c("term","gene")

# Prepare Msigdb geneset
MSIGDB_C3_TFT.GTRD <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%  dplyr::select(gs_name, entrez_gene)
MSIGDB_C3_TFT.TFT_Legacy <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>%  dplyr::select(gs_name, entrez_gene)
MSIGDB_C2_CP.BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>%  dplyr::select(gs_name, entrez_gene)
MSIGDB_H <- msigdbr(species = "Homo sapiens", category = "H") %>%  dplyr::select(gs_name, entrez_gene)
MSIGDB_C7_IMMUNESIGDB <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") %>%  
  dplyr::select(gs_name, entrez_gene) %>% 
  filter(grepl("NK", gs_name) ) #& !grepl("CD8|CD4|PDC|EPITHELIAL_CELLS|MACROPHAGE|MYELOID|MONO|NK|DC|ILC", gs_name))

#####################################################################################
## Prepare Background list for SYMBOL ID functionnal datasets ORA  ####
#####################################################################################
## @knitr format_background_data

# Background list conversion for ENTREID functional enrichment analysis (such enrichKEGG)
UNIVERSE_EntrezID <- AnnotationDbi::select(org.Hs.eg.db, UNIVERSE_SymbolID, 'ENTREZID','SYMBOL')
UNIVERSE_EntrezID <- UNIVERSE_EntrezID[which(UNIVERSE_EntrezID$ENTREZID != is.na(UNIVERSE_EntrezID$ENTREZID)),]
UNIVERSE_SymbolID <- UNIVERSE_EntrezID$SYMBOL
UNIVERSE_EntrezID <- UNIVERSE_EntrezID$ENTREZID

#######################
## Compare Cluster ####
#######################

# Here we perform cluster functional enrichment and compare them
## @knitr Compare_cluster_analysis

# Create an empty list to store the 9 ORA analysis
CompareCLUST <- vector(mode = "list", length = 9)
names(CompareCLUST)<-c("Staudt", "Reactome", "GO", "KEGG", "MSIGDB_TF_GTRD", "MSIGDB_TF_Legacy", "MSIGDB_BIOCARTA", "MSIGDB_HALLMARK","MSIGDB_IMMUNESIGDB")

# Compute the 9 ORA analysis
CompareCLUST[[1]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enrichGO", OrgDb = "org.Hs.eg.db", universe = UNIVERSE_EntrezID, readable = TRUE)
CompareCLUST[[2]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enrichPathway", universe = UNIVERSE_EntrezID,  readable=TRUE)
CompareCLUST[[3]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enrichGO", OrgDb = "org.Hs.eg.db",ont= "BP", universe = UNIVERSE_EntrezID, readable = TRUE)
CompareCLUST[[4]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enricher", TERM2GENE = MSIGDB_C3_TFT.GTRD, universe = UNIVERSE_EntrezID)
CompareCLUST[[5]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enricher", TERM2GENE = MSIGDB_C3_TFT.GTRD, universe = UNIVERSE_EntrezID)
CompareCLUST[[6]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enricher", TERM2GENE = MSIGDB_C3_TFT.TFT_Legacy, universe = UNIVERSE_EntrezID)
CompareCLUST[[7]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enricher", TERM2GENE = MSIGDB_C2_CP.BIOCARTA, universe = UNIVERSE_EntrezID)
CompareCLUST[[8]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enricher", TERM2GENE = MSIGDB_H, universe = UNIVERSE_EntrezID)
CompareCLUST[[9]] <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enricher", TERM2GENE = MSIGDB_C7_IMMUNESIGDB, universe = UNIVERSE_EntrezID)

# Set genes readable (SYMBOL)
for(j in 4:9){
  CompareCLUST[[j]] <- setReadable(CompareCLUST[[j]], OrgDb = org.Hs.eg.db, keyType="ENTREZID")
}


# NK vs NK2 vs NK3

#Perform GO analysis on a given seurat Object

#Part 1:
#  The aim of this script is to retrieve Functionnal Annotation Datasets
#  & Format both these datasets, DEG lists and Background Lists for ORA


library(clusterProfiler)
library(rlist)


PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")
Markers= FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = THR_log2FC , verbose = TRUE)
saveRDS(Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/DEG/DEG_NK_1ABC_2AB_3/MarkersNK1ABC_2AB_3.rds")

levels(PBMC$FirstClust) = c("NK1", "NK3", "NK1", "NK1", "NK2", "NK2")
PBMC = SetIdent(PBMC, value= "FirstClust")
DimPlot(PBMC)
Markers= FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = THR_log2FC , verbose = TRUE)
saveRDS(Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/DEG/DEG_NK1vsNK2vsNK3/MarkersNK1.2.3.rds")


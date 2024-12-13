#Reinject with all the genes
#OpenAll


PBMC_with_All_Genes = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/00_RawData/Seurat_Obj_other_samples/Sent_By_Janine/merged_and_allgenes.rds")

PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1.rds")


PBMC_with_All_Genes = subset(PBMC_with_All_Genes, cells= Cells(PBMC) )

PBMC_with_All_Genes@meta.data = PBMC@meta.data

saveRDS( PBMC_with_All_Genes , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")


PBMC_for_Metadata = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")
PBMC@meta.data = PBMC_for_Metadata@meta.data

PBMC@reductions = PBMC_for_Metadata@reductions

PBMC= SetIdent(PBMC, value= "FirstClust")

PBMC_Norm= as.SingleCellExperiment(PBMC)
PBMC_Norm = multiBatchNorm(PBMC_Norm,  batch = PBMC_Norm$orig.ident)
PBMC_Norm= as.Seurat(PBMC_Norm)

All_genes= rownames(PBMC_Norm)
PBMC=ScaleData(PBMC_Norm, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )

saveRDS(PBMC, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")


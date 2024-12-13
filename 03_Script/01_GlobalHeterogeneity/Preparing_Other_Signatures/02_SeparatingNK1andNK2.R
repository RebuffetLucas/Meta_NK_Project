#Separate Bright and DIm

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

#Merge NK1A and NK1B
PBMC= SetIdent(PBMC, value= "FirstClust")
table(PBMC$FirstClust)
levels(PBMC$FirstClust) = c("NK1B", "NK3",  "NK1A" ,"NK1A", "NK2B", "NK2A")
PBMC= SetIdent(PBMC, value= "FistClust")
table(PBMC$FirstClust)

PBMC_Keep = PBMC

PBMC_Bright = subset(PBMC, idents= c( "NK2A" , "NK2B"))
PBMC_Dim = subset(PBMC, idents= c("NK1A", "NK1B", "NK3"))


PBMC= PBMC_Bright
PBMC= PBMC_Dim
PBMC$FirstClust = droplevels(PBMC$FirstClust)
PBMC = SetIdent(PBMC, value= "FirstClust")
MONITORED_Markers= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Merging_NK1AandNK1B/TOP25_NK1A_Band3__5MainClusters.rds")

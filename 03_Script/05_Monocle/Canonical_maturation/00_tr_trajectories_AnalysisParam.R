# Enter the parameters of the analysis for tr_trajectories analysis

  
### For Dim ###
NK_Seurat=readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
NK_Seurat = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

palette<-c('NK1C'='#F8766D','NKint'='#8494FF',
           'NK1A'='#0CB702',
           'NK1B'='#00BFC4','NK2'='#ED68ED',
           'NK3'='#ABA300')

table(NK_Seurat$FirstClust)

NK_Seurat = SetIdent(NK_Seurat, value= "FirstClust")
DimPlot(NK_Seurat, cols= palette)


CLUSTERS_TO_REMOVE = c("NK3","NK2")

LIST_PARAM_CONTROL =list(ncenter=150)
root_pr_nodes = "Y_1"

#For graph test:

#Q_Value_Limit = 0.00005

#Coefficients
#k=9





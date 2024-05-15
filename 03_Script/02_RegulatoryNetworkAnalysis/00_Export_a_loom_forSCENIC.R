library(scater)
library(SeuratDisk)
library(loomR)
library(SCopeLoomR)
### This script creates a loom file for Dim , Bright and Dim+Bright+Prolif+Adapt
#It also allows to export the metadata table
#This version V2 works well for further analysis



#Create a loom for Bright:
NK_All=readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2Chem.rds")

DimPlot(NK_All, reduction = "UMAP")

build_loom(file.name = "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/data/NK_All/NKAll.loom",
           dgem= NK_All@assays$RNA@counts,
           title= "NK_All" ,
          default.embedding= NK_All@reductions$UMAP@cell.embeddings ,
           default.embedding.name="UMAP"
           )

loom <- open_loom("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/data/NK_All/NKAll.loom", mode = "r+")
close_loom(loom)

NK_All_info= NK_All@meta.data

saveRDS(NK_All_info , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/02_RegulatoryNetworkAnalysis/data/NK_All/NK_All_Info.rds")











#Create a loom for All_NK (all genes)
NK_All=readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/01_GlobalHeterogeneity/Output_First_Clustering/ChemData100/ChemData_Rescaled_Norm100/With_NK3_Complete/Rename_imm/allNKallgenes.rds")



NK_All= SetIdent(NK_All, value = "SecondClust" ) 
NK_BrightDimAdaptProlif= subset(NK_All, idents=c("NK_Dying"), invert= TRUE )
NK_BrightDimAdaptProlif$SecondClust = droplevels(NK_BrightDimAdaptProlif$SecondClust)
NK_BrightDimAdaptProlif= SetIdent(NK_BrightDimAdaptProlif, value = "SecondClust" ) 
  
NK_All@reductions$UMAP@cell.embeddings = NK_All2@reductions$UMAP@cell.embeddings

DimPlot(NK_All, reduction = "umap")
DimPlot(NK_All, reduction = "UMAP", group.by = "SecondClust")

NK_All2@meta.data= NK_All@meta.data

dim(NK_All2)
build_loom(file.name = "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom",
           dgem= NK_BrightDimAdaptProlif@assays$RNA@counts,
           title= "NK_All"
)

NK_All_info= NK_BrightDimAdaptProlif@meta.data
saveRDS(NK_All_info , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK_Info.rds")

loom <- open_loom("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom", mode = "r+")
close_loom(loom)



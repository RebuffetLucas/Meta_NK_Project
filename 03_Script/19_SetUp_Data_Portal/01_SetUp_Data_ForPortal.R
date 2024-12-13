
library(SeuratWrappers)
# devtools::install_github("hhoeflin/hdf5r")

PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")
PBMC = SetIdent(PBMC, value= "FirstClust")


#Remove stuff
  #Remove reductions
PBMC@reductions[["HARMONY"]] = NULL
PBMC@reductions[["TSNE"]] = NULL
PBMC@reductions[["PCA"]] =  NULL

  #Remove command
PBMC@commands[["ScaleData.RNA"]] =  NULL

colnames(PBMC@meta.data)

PBMC@meta.data[["Chemistry"]] = NULL
PBMC@meta.data[["nCount_RNA"]] = NULL
PBMC@meta.data[["nFeature_RNA"]] = NULL
PBMC@meta.data[["percent.mito"]] = NULL
PBMC@meta.data[["percent.ribo"]] = NULL
PBMC@meta.data[["orig.ident"]] = NULL
PBMC@meta.data[["Dataset"]] = NULL
PBMC@meta.data[["donor"]] = NULL
PBMC@meta.data[["sizeFactor"]] = NULL
PBMC@meta.data[["code"]] = NULL
PBMC@meta.data[["cluster"]] = NULL
PBMC@meta.data[["FirstClust"]] = NULL
PBMC@meta.data[["nkg2c"]] = NULL
PBMC@meta.data[["SecondClust"]] = NULL
PBMC@meta.data[["Crinier_13_hNKGenes_20181"]] = NULL


#Have a look
FEATURE= "NEAT1"
GROUP = "FirstClust"

FeaturePlot(PBMC, features= FEATURE)
DimPlot(PBMC, group.by= GROUP)

VlnPlot(PBMC, features= FEATURE)
VlnPlot(PBMC, features= FEATURE, pt.size = 0)

PBMC$ident = as.character(PBMC$ident)


DimPlot(PBMC,  cols= palette) & NoAxes()

table(PBMC$ident)

#Convert to AnnData and Save withtout regulons

SaveH5Seurat(PBMC, filename = "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_Noreg.h5Seurat", overwrite = TRUE)
Convert("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_Noreg.h5Seurat", dest = "h5ad", overwrite = TRUE)


#Add the Regulon slot

regulonAUC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/regulonAUC.rds")

cells = colnames(regulonAUC)


PBMC2 = subset(PBMC, cells = cells) #Keep only cells that have regulons


#Save without Regulons but only with less cells that have regulons
SaveH5Seurat(PBMC2, filename = "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_reduced_Noreg.h5Seurat")
Convert("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_reduced_Noreg.h5Seurat", dest = "h5ad")



PBMC2


To_Insert = regulonAUC@assays@data@listData[["AUC"]]
rownames(To_Insert) = gsub('.{3}$', '', rownames(To_Insert))

#Create Regulon assay
Regulons_Insert = CreateAssayObject(counts = To_Insert)

#Create regulon slot in seurat object
PBMC2@assays[["Regulons"]] = Regulons_Insert

#test that it's working well
sum(is.na(PBMC2@assays[["Regulons"]]@counts))
sum(is.finite(PBMC2@assays[["Regulons"]]@counts))
sum(!is.finite(PBMC2@assays[["Regulons"]]@counts))
all(PBMC2@assays[["Regulons"]]@counts == PBMC2@assays[["Regulons"]]@data)


#PBMC2@assays[["Regulons"]]@data = PBMC2@assays[["Regulons"]]@counts  # Not sure this is usefull

#Have a look at the slot
PBMC2@assays[["RNA"]]@data

#Check the slot can be used
DefaultAssay(PBMC2) = "Regulons"
FeaturePlot(PBMC2, features = "TCF7")

PBMC2@reductions[["UMAP"]]@misc  = list() #This fix does not work






#Convert to AnnData and Save with regulons

saveRDS(PBMC2 , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_WithReg.rds" )
SaveH5Seurat(PBMC2, filename = "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_WithReg.h5Seurat")
Convert("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_WithReg.h5Seurat", dest = "h5ad")

write.csv(regulonAUC@assays@data@listData[["AUC"]] , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/regulonAUC.csv")

dim(regulonAUC@assays@data@listData[["AUC"]])
test_csv = read_csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/regulonAUC.csv")


#Save one with Regulons only
PBMC2@assays[["RNA"]] =  Regulons_Insert

SaveH5Seurat(PBMC2, filename = "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_RegOnly.h5Seurat", overwrite = TRUE)
Convert("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/19_SetUp_Data_Portal/PBMC_WithReg.h5Seurat", dest = "h5ad")






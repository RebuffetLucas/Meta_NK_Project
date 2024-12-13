#Quick look at Tang data

library(SeuratDisk)


#Blood Only
Convert("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5ad", dest="h5seurat", overwrite = TRUE) # Creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
PBMC_Blood = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5seurat")

#All derived
Convert("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5ad", dest="h5seurat", overwrite = TRUE) # Creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
PBMC_All = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5seurat")



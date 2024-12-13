#Removing the cells that are not interesting

#Try to cluster them only in 3 pops then look at their signatures

Seurat_NK = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")
#Here we choose not to remove D3 and D7 and we will precise D0 as reference during the integration process

#Remove previous clasification to make sure
Seurat_NK@reductions[["apca"]] = NULL
Seurat_NK@reductions[["aumap"]] = NULL
Seurat_NK@reductions[["pca"]] = NULL
Seurat_NK@reductions[["spca"]] = NULL
Seurat_NK@reductions[["umap"]] = NULL
Seurat_NK@reductions[["wnn.umap"]] = NULL



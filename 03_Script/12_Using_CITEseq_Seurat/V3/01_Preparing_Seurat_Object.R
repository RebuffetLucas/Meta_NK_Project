#Removing the cells that are not interesting

#Try to cluster them only in 3 pops then look at their signatures

Seurat_NK = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")
#Seurat_NK = subset(Seurat_NK, idents= "NK Proliferating" , invert= TRUE)
Seurat_NK = subset(Seurat_NK, subset= time =="0")

#Remove previous clasification to make sure
Seurat_NK@reductions[["apca"]] = NULL
Seurat_NK@reductions[["aumap"]] = NULL
Seurat_NK@reductions[["pca"]] = NULL
Seurat_NK@reductions[["spca"]] = NULL
Seurat_NK@reductions[["umap"]] = NULL
Seurat_NK@reductions[["wnn.umap"]] = NULL

Seurat_NK@assays[["SCT"]]@SCTModel.list= NULL


#Seurat_NK2 = UpdateSeuratObject(Seurat_NK)
#Seurat_NK2 = UpdateSCTAssays(Seurat_NK2)



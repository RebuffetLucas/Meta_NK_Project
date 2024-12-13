#Dependencies

library(SeuratDisk)

reference <- LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/pbmc_multimodal.h5seurat")

reference = SetIdent(reference,  value= "celltype.l1")

table(reference$celltype.l1)
table(reference$celltype.l2)
table(reference$celltype.l3)


table(reference$donor)
table(reference$time)
table(reference$lane)

Seurat_NK = subset(reference, idents= "NK")

DimPlot(Seurat_NK, group.by= "celltype.l3", reduction = "wnn.umap")

Seurat_NK = SetIdent(Seurat_NK,  value= "celltype.l3")
table(Seurat_NK$celltype.l3)

Seurat_NK$celltype.l3 = droplevels(Seurat_NK$celltype.l3)

saveRDS(Seurat_NK, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")



Seurat_NK_Final = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_After_Satija_Procedure.rds")

Seurat_NK_Final = SetIdent(Seurat_NK_Final , value = "wsnn_res.0.3")

DimPlot(Seurat_NK_Final, group.by = "wsnn_res.0.6" )

DimPlot(Seurat_NK_Final)

table(Seurat_NK_Final$wsnn_res.0.3)

Seurat_NK_Final = subset(Seurat_NK_Final, idents= c("0","1","2","3","4"))

Seurat_NK_Final <- FindMultiModalNeighbors(
  Seurat_NK_Final, reduction.list = list("r_pca", "r_apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)


#### Downstream analysis ####

Seurat_NK_Final  = RunUMAP(Seurat_NK_Final, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_NK_Final <- FindClusters(Seurat_NK_Final, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = TRUE)



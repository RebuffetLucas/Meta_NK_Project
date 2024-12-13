
#Basic treatment on Seurat_NK2

#Initialize the object
Seurat_NK2= SCTransform(Seurat_NK, assay= "SCT", new.assay.name = "SCT2" ) #This step includes Normalize + scale + find variable features

Seurat_NK2@assays[["SCT"]] = NULL
DefaultAssay(Seurat_NK2) <- 'SCT2'


Seurat_NK2 <- RunPCA(Seurat_NK2, assay = "SCT2",verbose = FALSE)


#Seurat_NK2 <- NormalizeData(Seurat_NK2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(Seurat_NK2) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(Seurat_NK2) <- rownames(Seurat_NK2[["ADT"]])
Seurat_NK2 <- NormalizeData(Seurat_NK2, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

Seurat_NK2 <- FindMultiModalNeighbors(
  Seurat_NK2, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

Seurat_NK2 <- RunUMAP(Seurat_NK2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_NK2 <- FindClusters(Seurat_NK2, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = FALSE)

p1 <- DimPlot(Seurat_NK2, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) 
p2 <- DimPlot(Seurat_NK2, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) 
p3 <- DimPlot(Seurat_NK2, reduction = 'wnn.umap', group.by = 'celltype.l3', label = TRUE, repel = TRUE, label.size = 2.5) 
p4 <- DimPlot(Seurat_NK2, reduction = 'wnn.umap', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5)
p1 + p2 +p3 + p4

#Look at batch effects
p1 <- DimPlot(Seurat_NK2, reduction = 'apca', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5)
p2  <- DimPlot(Seurat_NK2, reduction = 'pca', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5)
p1 + p2

p1 <- DimPlot(Seurat_NK, reduction = 'apca', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5)
p2  <- DimPlot(Seurat_NK, reduction = 'pca', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5)
p1 + p2



table(Seurat_NK2$lane)

DimPlot(Seurat_NK2, group.by = "time", reduction = "wnn.umap")  

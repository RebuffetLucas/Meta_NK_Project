

#### Dealing with RNA first ####
Seurat_NK2= SCTransform(Seurat_NK, assay= "SCT", new.assay.name = "SCT" )
#Seurat_NK2@assays[["ADT"]] = NULL
DefaultAssay(Seurat_NK2) <- 'SCT'
#Partition by patients and time point (here we only have D0)
Seurat_NK_List = SplitObject(Seurat_NK2, split.by = "orig.ident")
#SCT on the separated elements
Seurat_NK_List = lapply(X= Seurat_NK_List, FUN= function(x){
  x = SCTransform(x, assay= "SCT", new.assay.name = "SCT", do.scale = DATA_SCALE, do.center = DATA_CENTER, return.only.var.genes = FALSE )
})

#Apply rPCA
#Integration features
variable_RNA_features = SelectIntegrationFeatures(object.list = Seurat_NK_List, nfeatures = 3000 )
Seurat_NK_List <- PrepSCTIntegration(object.list = Seurat_NK_List, anchor.features = variable_RNA_features)

#Run PCA
Seurat_NK_List = lapply(X= Seurat_NK_List, FUN= function(x){
  x = RunPCA(x, assay= "SCT" , features= variable_RNA_features )
})

#Find anchors
immune.anchors <- FindIntegrationAnchors(object.list = Seurat_NK_List, normalization.method = "SCT",
                                         anchor.features = variable_RNA_features, dims = 1:30, reduction = "rpca", k.anchor = 20)
#Create the integrated object
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE, reduction.name= "r_pca")

#Have a look
p1 = DimPlot(immune.combined.sct, group.by = "orig.ident")


DimPlot(immune.combined.sct, group.by = "orig.ident")



#### Dealing with ADT ####
DefaultAssay(immune.combined.sct) <- 'ADT'

#All ADT are considered as variable
VariableFeatures(immune.combined.sct) <- rownames(immune.combined.sct[["ADT"]])

#Normalize and scale  
immune.combined.sct <- NormalizeData(immune.combined.sct, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData( do.scale = DATA_SCALE, do.center = DATA_CENTER)

#Split for rPCA on ADT
Seurat_NK_List = SplitObject(immune.combined.sct, split.by = "orig.ident")

#All ADT are considered as variable
variable_ADT_features = VariableFeatures(immune.combined.sct) 

#Run PCA
Seurat_NK_List = lapply(X= Seurat_NK_List, FUN= function(x){
  x = RunPCA(x, assay= "ADT" , features= variable_ADT_features, reduction.name= "pca" )
})

#Find anchors
immune.anchors2 <- FindIntegrationAnchors(object.list = Seurat_NK_List,   assay= NULL, normalization.method = "LogNormalize",
                                          anchor.features = variable_ADT_features, dims = 1:30, reduction = "rpca", k.anchor = 20)

immune.combined.sct.adt <- IntegrateData(anchorset = immune.anchors2, normalization.method = "LogNormalize", dims = 1:30, new.assay.name = "integratedADT")
immune.combined.sct.adt <- ScaleData(immune.combined.sct.adt, do.scale = DATA_SCALE, do.center = DATA_CENTER, verbose = TRUE)
immune.combined.sct.adt <- RunPCA(immune.combined.sct.adt, verbose = FALSE, reduction.name= "r_apca" )

p2 = DimPlot(immune.combined.sct.adt , reduction = "r_apca", group.by = "orig.ident")
DimPlot(immune.combined.sct.adt , reduction = "r_apca", group.by = "orig.ident")

p1 + p2


#Go back to a simple Seurat object with my r-pca and my r_apca

Seurat_NK_Final = immune.combined.sct.adt
Seurat_NK_Final@reductions[["r_pca"]] = immune.combined.sct[["r_pca"]]


#### Run WNN on rPCA and rAPCA ####

DefaultAssay(Seurat_NK_Final) <- 'ADT'


Seurat_NK_Final <- FindMultiModalNeighbors(
  Seurat_NK_Final, reduction.list = list("r_pca", "r_apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)


#### Downstream analysis ####

Seurat_NK_Final  = RunUMAP(Seurat_NK_Final, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_NK_Final <- FindClusters(Seurat_NK_Final, graph.name = "wsnn", algorithm = 3, resolution = 0.4, verbose = TRUE)

DimPlot(Seurat_NK_Final)
#Saving
saveRDS(Seurat_NK_Final, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_After_Satija_Procedure_OnlyD0.rds")

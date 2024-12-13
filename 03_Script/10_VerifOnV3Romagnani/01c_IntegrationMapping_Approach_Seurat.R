#Use V2 as a reference

PBMC_V2 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
DimPlot(PBMC_V2, label = TRUE)

#Load V3 samples

#PBMC_V3 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/00_RawData/Seurat_Obj_other_samples/Sent_By_Janine/merged_and_downsampled.rds")
#PBMC_V3 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/00_RawData/seurat_after_manual_filteringAndchangeMetadata/IPH_NK5.rds")

#PBMC_V3= subset(PBMC_V3, subset= Chemistry == "V3" )
#PBMC_V3 = subset(PBMC_V3, features= rownames(PBMC_V2))

All_genes= rownames(PBMC_V3)
PBMC_V3=ScaleData(PBMC_V3, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )


#Find the variable features of V2 for transfer
#Merge them
Merged_Seurat= PBMC_V2
Merged_Seurat$orig.ident = droplevels(Merged_Seurat$orig.ident)


#Free Memory
gc()


#Rescale with MultiBatchNorm (take into account the sequencing depth of each Dataset)
# !!! This step replace NormalizeDataStep !!!
sce_Merged= as.SingleCellExperiment(Merged_Seurat)
sce_Merged_Rescaled = multiBatchNorm(sce_Merged,  batch = sce_Merged$orig.ident)
Merged_Seurat_Rescaled= as.Seurat(sce_Merged_Rescaled)

#Only LogNorm
#Merged_Seurat_Rescaled = NormalizeData(Merged_Seurat) #Do not run this step if MultiBatchNorm was used

#Find VariableFeatures across samples
#Variable on orig.ident
seurat_list_Rescaled = SplitObject(Merged_Seurat_Rescaled, split.by = "orig.ident")
seurat_list_Rescaled <- lapply(seurat_list_Rescaled, function(x) FindVariableFeatures(x, nfeatures= 5000, selection.method = "vst"))
VARIABLE_FEATURES = SelectIntegrationFeatures(  seurat_list_Rescaled,  nfeatures = 2000,  assay = NULL,  verbose = TRUE)
#Free Memory
rm(sce_Merged, sce_Merged_Rescaled, seurat_list_Rescaled)
gc()


#Transfer
PBMC.anchors <- FindTransferAnchors(reference = PBMC_V2, query = PBMC_V3,
                                        dims = 1:30, reference.reduction = "pca", features = VARIABLE_FEATURES)
predictions <- TransferData(anchorset = PBMC.anchors, refdata = PBMC_V2$FirstClust,
                            dims = 1:30)

PBMC_V3 <- AddMetaData(PBMC_V3, metadata = predictions)

#Have a look at the results
table(PBMC_V3$predicted.id)

VlnPlot(PBMC_V3 , c("ACTB","SPON2",  "ACTG1", "PRF1", "PFN1"), group.by = "predicted.id")
VlnPlot(PBMC_V3 , c("KLRC2","GZMH",  "IL32", "CD52"), group.by = "predicted.id", pt.size = 0)

VizDimLoadings(PBMC_V2 , reduction = "pca" , dims= 1:2)
VizDimLoadings(PBMC_V3 , reduction = "pca" , dims= 1:2)



#Unimodal UMAP Projection
PBMC_V2 <- RunUMAP(PBMC_V2, dims = 1:30, reduction = "harmony", return.model = TRUE)
PBMC_V3 <- MapQuery(anchorset = PBMC.anchors, reference = PBMC_V2, query = PBMC_V3, refdata = list(FirstClust = "FirstClust"), reference.reduction = "harmony", reduction.model = "umap")

p1 <- DimPlot(PBMC_V2, reduction = "umap", group.by = "FirstClust", label = TRUE, label.size = 3,
              repel = TRUE, cols = palette) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(PBMC_V3, reduction = "ref.umap", group.by = "predicted.FirstClust", label = TRUE,
              label.size = 3, repel = TRUE, cols = palette, pt.size = 0.6) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(PBMC_V3, reduction = "ref.umap", group.by = "nkg2c", label = F,
              label.size = 3, repel = TRUE ,pt.size = 0.6) + ggtitle("Query transferred labels")

#Look at differencialy expressed genes
PBMC_V3 = SetIdent(PBMC_V3, value= "predicted.FirstClust")

All_Markers = FindAllMarkers(PBMC_V3 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


#Look at markers like Constance
markers_overCl=All_Markers

topGenes=c()
for (cl in sort(unique(Idents(PBMC_V3)))) {
  print(cl)
  markers_overCl_cl=markers_overCl[
    markers_overCl$avg_log2FC>0&
      markers_overCl$p_val_adj<=0.05&
      markers_overCl$cluster==cl,
  ]
  
  markers_overCl_cl=	
    markers_overCl_cl[order(markers_overCl_cl$p_val),]
  topGenes_cl=markers_overCl_cl$gene[
    1:min(FINDMARKERS_SHOWTOP,nrow(markers_overCl_cl))]
  topGenes[[cl]]=topGenes_cl
  print(	topGenes_cl)
}

topGenes2=topGenes

topGenesb=unique(unlist(topGenes))

DotPlot(PBMC_V3, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

DotPlot(PBMC_V3, features = topGenesb , cols = "Spectral", col.max = 1.5 , col.min = -1.5) + theme(axis.text.x = element_text(angle = 90))




DimPlot(PBMC_V2 , reduction = "umap")



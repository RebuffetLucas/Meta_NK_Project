
#Load Ref and Normalize ADT
reference = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")

reference= subset(reference, idents= "NK Proliferating", invert=TRUE)
DefaultAssay(object = reference) <- "ADT"
reference = NormalizeData(reference, normalization.method = 'CLR', margin = 2)

DefaultAssay(object = reference) <- "SCT"
DimPlot(reference, reduction = "wnn.umap")

#LoadQuery
query = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK4/05_Output/01_GlobalHeterogeneity/PBMC_All.rds")

query_Metadata = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
query= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")

query@meta.data = query_Metadata@meta.data

query =  SCTransform(query, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)


query = MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l2",
    celltype.l3 = "celltype.l3",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 = DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)  
p2 = DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE)  
p1 + p2


DimPlot(query, reduction = "ref.umap", group.by = "FirstClust", label = TRUE, label.size = 3, repel = TRUE)  
DimPlot(reference, reduction = "spca", label = TRUE, label.size = 3, repel = TRUE)  
DimPlot(reference, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE)  


#Have a look at the predicted ADT

table(query$predicted.celltype.l2, query$FirstClust)


DefaultAssay(query) = "predicted_ADT"
query = NormalizeData(query, normalization.method = 'CLR', margin = 2)

query = SetIdent(query, value= "FirstClust")

Markers_Infered_ADT =  FindAllMarkers(query , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

Markers_Infered_ADT =  FindAllMarkers(query , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , logfc.threshold = FINDMARKERS_LOGFC_THR, verbose = TRUE)

query@reductions[["UMAP"]]=query_Metadata@reductions[["UMAP"]]


FeaturePlot(query, reduction = "UMAP", features= Markers_Infered_ADT$gene[1:9] )

p1 = FeaturePlot(query, reduction= "UMAP", features= unique(Markers_Infered_ADT$gene)[1:9], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p2 = FeaturePlot(query, reduction= "UMAP", features= unique(Markers_Infered_ADT$gene)[10:18], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p3 = FeaturePlot(query, reduction= "UMAP", features= unique(Markers_Infered_ADT$gene)[18:24], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

p1
p2
p3


p5 = VlnPlot(query,   features= unique(Markers_Infered_ADT$gene)[1:9], pt.size=0)
p6 = VlnPlot(query,   features= unique(Markers_Infered_ADT$gene)[10:18], pt.size=0) 
p7 = VlnPlot(query,   features= unique(Markers_Infered_ADT$gene)[18:24], pt.size=0)

p9 = VlnPlot(query,   features= c("CD161"), pt.size=0)

p10 = VlnPlot(query,   features= c("CD158"), pt.size=0)

p9

p5
p6
p7
p8



#DotPlot(query, features= Markers_Infered_ADT$gene, pt.size = 0)
DotPlot(query, features = unique(Markers_Infered_ADT$gene) , cols = "Spectral", col.max = 1 , col.min = 0) + theme(axis.text.x = element_text(angle = 90))
DotPlot(query, features = unique(Markers_Infered_ADT$gene) , cols = "Spectral", col.max = 1 , col.min = 0) + theme(axis.text.x = element_text(angle = 90))
DotPlot(query, features = unique(Markers_Infered_ADT$gene) , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))
DotPlot(query, features = unique(Markers_Infered_ADT$gene) , cols = "Spectral", scale=FALSE) + theme(axis.text.x = element_text(angle = 90))

VlnPlot(query, features = unique(Markers_Infered_ADT$gene)[1:9], pt.size=0)

DimPlot(query, reduction= "UMAP", group.by = "FirstClust", label = TRUE)



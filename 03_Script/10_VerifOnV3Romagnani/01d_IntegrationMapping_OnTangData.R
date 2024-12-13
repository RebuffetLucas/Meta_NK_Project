#Importing Tang et al data

#Cells_Tang= Read10X("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/Zhang_Review_Cell/row_names/", gene.column = 2, cell.column = 2)
#Cells_Tang= Read10X("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/Zhang_Review_Cell/row_names/")


mat= readMM("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/Zhang_Review_Cell/row_names2/NK_rawcount.mtx")
genes= read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/Zhang_Review_Cell/row_names2/NK_rawcount.genes")
barcodes=read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/Zhang_Review_Cell/row_names2/NK_rawcount.barcodes")

dim(matrix)
dim(genes)
dim(barcodes)

colnames(mat)=genes[,2]
rownames(mat)=barcodes[,2]
Tang_Object = CreateSeuratObject(t(mat) , min.cells = 3 , min.features = 200)
Tang_Metadata = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/Zhang_Review_Cell/NK_rawcount_metadata.csv")

rownames(Tang_Metadata) = Tang_Metadata$X

Tang_Object2 = AddMetaData(Tang_Object, metadata =Tang_Metadata )

#Which metadata to regress on?
VlnPlot(Tang_Object2, features = c("nCount_RNA", "nFeature_RNA"), group.by = "sampleID", pt.size = 0)
VlnPlot(Tang_Object2, features = c("nCount_RNA", "nFeature_RNA"), group.by = "meta_patientID", pt.size = 0)
VlnPlot(Tang_Object2, features = c("nCount_RNA", "nFeature_RNA"), group.by = "meta_histology", pt.size = 0)
VlnPlot(Tang_Object2, features = c("nCount_RNA", "nFeature_RNA"), group.by = "meta_histology", pt.size = 0)
VlnPlot(Tang_Object2, features = c("nCount_RNA", "nFeature_RNA"), group.by = "meta_tissue", pt.size = 0)

#Reference Dataset preprocessing

PBMC_V2 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
DimPlot(PBMC_V2, label = TRUE)


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

PBMC_V2 <- RunUMAP(PBMC_V2, dims = 1:30, reduction = "harmony", return.model = TRUE)



#QueryDataset_Preprocessing
table(Tang_Object2$meta_histology)
Tang_Object2 = subset(Tang_Object2, features= rownames(PBMC_V2))
Query_Dataset = Tang_Object2


PBMC_Norm= as.SingleCellExperiment(Query_Dataset)
PBMC_Norm = multiBatchNorm(PBMC_Norm,  batch = PBMC_Norm$meta_histology) #ADJUST_PARAM_HERE
Query_Dataset= as.Seurat(PBMC_Norm)

All_genes= rownames(Query_Dataset)
Query_Dataset=ScaleData(Query_Dataset, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )
Query_Dataset_batches = SplitObject(Query_Dataset, split.by = "meta_histology") #ADJUST_PARAM_HERE





#No scaling
#All_genes= rownames(Query_Dataset)
#Query_Dataset=ScaleData(Query_Dataset, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER ) 


#Mapping
  #Individual anchors
anchors <- list()
for (i in 1:length(Query_Dataset_batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = PBMC_V2,
    query = Query_Dataset_batches[[i]],
    k.filter = NA,
    reference.reduction = "pca", 
    #reference.neighbors = "FirstClust", 
    dims = 1:30,
    features = VARIABLE_FEATURES
  )
}

  #Individual Mapping
for (i in 1:length(Query_Dataset_batches)) {
  Query_Dataset_batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = Query_Dataset_batches[[i]],
    reference = PBMC_V2, 
    refdata = list(
      FirstClust = "FirstClust" 
      ),
    reference.reduction = "harmony",
    reduction.model = "umap"
  )
}

#Don't do the cash step but can be done to accelerate:
# bm <- FindNeighbors(
#object = bm,
#reduction = "spca",
#dims = 1:50,
#graph.name = "spca.annoy.neighbors", 
#k.param = 50,
#cache.index = TRUE,
#return.neighbor = TRUE,
#l2.norm = TRUE )



#Explore the mapping results


p1 <- DimPlot(Query_Dataset_batches[[1]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[1]]$meta_histology)))
p2 <- DimPlot(Query_Dataset_batches[[2]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[2]]$meta_histology)))
p3 <- DimPlot(Query_Dataset_batches[[3]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[3]]$meta_histology)))
p4 <- DimPlot(Query_Dataset_batches[[4]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[4]]$meta_histology)))
p5 <- DimPlot(Query_Dataset_batches[[5]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[5]]$meta_histology)))
p6 <- DimPlot(Query_Dataset_batches[[6]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[6]]$meta_histology)))
p7 <- DimPlot(Query_Dataset_batches[[7]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[7]]$meta_histology)))
p8 <- DimPlot(Query_Dataset_batches[[8]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette) + ggtitle(as.character(unique(Query_Dataset_batches[[8]]$meta_histology)))

p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8  + plot_layout(guides = "collect")


DimPlot(Query_Dataset_batches[[1]], reduction = 'ref.umap', group.by = 'predicted.FirstClust', label.size = 3, cols = palette, pt.size = 0.8) + ggtitle(as.character(unique(Query_Dataset_batches[[1]]$meta_histology)))

DimPlot(Query_Dataset_batches[[1]], reduction = 'ref.umap', group.by = 'Majortype', label.size = 3,  pt.size = 0.8) 
DimPlot(Query_Dataset_batches[[1]], reduction = 'ref.umap', group.by = 'celltype', label.size = 3,  pt.size = 0.8) 


#Merge the batches

Query_Dataset2 =merge(Query_Dataset_batches[[1]], Query_Dataset_batches[2:length(Query_Dataset_batches)] ,  merge.dr = "ref.umap") 
DimPlot(Query_Dataset2, reduction = "ref.umap", group.by =  "predicted.FirstClust", label = TRUE, repel = TRUE, label.size = 3, cols = palette) + NoLegend()


#See what is shared
Query_Dataset2$seurat_clusters = Query_Dataset2$predicted.FirstClust
p5 = ggplot(Query_Dataset2@meta.data, aes(x=Majortype, fill= seurat_clusters )) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Query_Dataset2@meta.data, aes(x=celltype, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p5 + p6

Query_Dataset2$seurat_clusters = Query_Dataset2$predicted.FirstClust
p5 = ggplot(Query_Dataset2@meta.data, aes(x=seurat_clusters, fill= Majortype )) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Query_Dataset2@meta.data, aes(x=seurat_clusters, fill=celltype)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p5 + p6

VlnPlot(Query_Dataset2, features = "predicted.FirstClust.score")

#Keeping only the reliable ones
Query_Dataset3 = subset(Query_Dataset2, predicted.FirstClust.score >0.6)

Query_Dataset3$seurat_clusters = Query_Dataset3$predicted.FirstClust
p5 = ggplot(Query_Dataset3@meta.data, aes(x=Majortype, fill= seurat_clusters )) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Query_Dataset3@meta.data, aes(x=celltype, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p5 + p6



Query_Dataset3$seurat_clusters = Query_Dataset3$predicted.FirstClust
p5 = ggplot(Query_Dataset3@meta.data, aes(x=seurat_clusters, fill= Majortype )) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Query_Dataset3@meta.data, aes(x=seurat_clusters, fill=celltype)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p5 + p6

DimPlot(Query_Dataset3, reduction = "ref.umap", group.by =  "predicted.FirstClust", label = F, repel = TRUE, label.size = 3, cols = palette, pt.size = 1) 



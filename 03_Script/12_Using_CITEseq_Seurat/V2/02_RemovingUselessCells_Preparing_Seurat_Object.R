#Removing the cells that are not interesting

#Try to cluster them only in 3 pops then look at their signatures

Seurat_NK = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/Seurat_CITEseq_NK.rds")
#Seurat_NK = subset(Seurat_NK, idents= "NK Proliferating" , invert= TRUE)
#Seurat_NK = subset(Seurat_NK, subset= time =="0")

#Remove previous clasification to make sure
Seurat_NK@reductions[["apca"]] = NULL
Seurat_NK@reductions[["aumap"]] = NULL
Seurat_NK@reductions[["pca"]] = NULL
Seurat_NK@reductions[["spca"]] = NULL
Seurat_NK@reductions[["umap"]] = NULL
Seurat_NK@reductions[["wnn.umap"]] = NULL

Seurat_NK@assays[["SCT"]]@SCTModel.list= NULL


Seurat_NK2 = UpdateSeuratObject(Seurat_NK)
Seurat_NK2 = UpdateSCTAssays(Seurat_NK2)

Seurat_NK2= SCTransform(Seurat_NK, assay= "SCT", new.assay.name = "SCT" )

#Rerun scaling and shits
DefaultAssay(Seurat_NK) <- 'SCT'
Seurat_NK <- NormalizeData(Seurat_NK) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(Seurat_NK) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(Seurat_NK) <- rownames(Seurat_NK[["ADT"]])
Seurat_NK <- NormalizeData(Seurat_NK, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

Seurat_NK <- FindMultiModalNeighbors(
  Seurat_NK, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

Seurat_NK <- RunUMAP(Seurat_NK, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_NK <- FindClusters(Seurat_NK, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = FALSE)

p1 <- DimPlot(Seurat_NK, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(Seurat_NK, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2


p1 <- DimPlot(Seurat_NK, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(Seurat_NK, reduction = 'wnn.umap', group.by = 'orig.ident', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

p1 <- DimPlot(Seurat_NK, reduction = 'wnn.umap', group.by = 'donor' , label = TRUE, repel = TRUE, label.size = 2.5)  
p2 <- DimPlot(Seurat_NK, reduction = 'wnn.umap', group.by = 'time', label = TRUE, repel = TRUE, label.size = 2.5)  
p1 + p2



df= as.data.frame(table(Idents(Seurat_NK)))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

p4 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    


p5 = ggplot(Seurat_NK@meta.data, aes(x=orig.ident, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p6 = ggplot(Seurat_NK@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p7 = ggplot(Seurat_NK@meta.data, aes(x=donor, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  

p4 + p5 + p6  + p7

p4 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    


p5 = ggplot(Seurat_NK@meta.data, aes(x=time, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p7 = ggplot(Seurat_NK@meta.data, aes(x=donor, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  

p5   + p7



table(Seurat_NK$lane)

DimPlot(Seurat_NK, group.by = "time", reduction = "wnn.umap")  

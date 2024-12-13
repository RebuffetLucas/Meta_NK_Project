#Scoring based on stress induced by Enzyme treatment

genes_of_interest= read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/coregene_df-FALSE-v3.csv")

genes_pos = genes_of_interest %>% filter(logFC>0)
list_genes= genes_pos$gene_symbol
list_genes= list_genes[!is.na(list_genes)]

list_genes = intersect(rownames(PBMC), list_genes)


Merged_Seurat_Rescaled = AddModuleScore(Merged_Seurat_Rescaled, features = list(list_genes) , name= "stress_score" , seed=19)

Merged_Seurat_Rescaled@meta.data <- Merged_Seurat_Rescaled@meta.data %>% dplyr::select(-matches("^Cluster"))


VlnPlot(Merged_Seurat_Rescaled, feature= "stress_score1", pt.size = 0, group.by = "Dataset")
VlnPlot(Merged_Seurat_Rescaled, feature= "stress_score1", pt.size = 0)

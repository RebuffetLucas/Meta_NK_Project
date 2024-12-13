#Create A umap OF THE USUAL SUSPECTS

library(rlist)
library(reshape2)




#Load Seurat Object
PBMC_for_Metadata = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")

PBMC@meta.data = PBMC_for_Metadata@meta.data

PBMC_Norm= as.SingleCellExperiment(PBMC)
PBMC_Norm = multiBatchNorm(PBMC_Norm,  batch = PBMC_Norm$orig.ident)
PBMC_Norm= as.Seurat(PBMC_Norm)


All_genes= rownames(PBMC_Norm)
PBMC=ScaleData(PBMC_Norm, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )



#Get the list of genes present in at list X percent of the cells


MIN_PCT_GENE= 5

#NK2

PBMC= SetIdent(PBMC, value = "FirstClust")

PBMC$FirstClust= droplevels(PBMC$FirstClust)



List_BrightvsDim = c("XCL1", "XCL2", "GZMK" , "FGFBP2", "GZMB", "KLRC1", "FCGR3A")



List_Genes1 = c(  "KLF2", "GZMA", "LGALS1", "ZEB2", "SELL", "CMA1", "ITGA4", "PRF1", "IL18RAP", "EOMES","SCIMP", "ITGAM", "ITGA2", "PLEK", "CD247" )
List_Genes2 = c("ITM2C", "RBPJ", "TNFSF10", "CD226", "TMEM176A", "TMEM176B", "IL7R","S1PR5", "CXCR6", "SDC4", "RGS1")



List_Markers = list.append(List_Genes1, List_Genes2)
List_of_Lists= list(List_Genes1, List_Genes2)

List_total= intersect(List_Markers, rownames(PBMC))


#Keep Genes expressed in more than X percent cells
perc_exp <- DotPlot(PBMC, features=List_total, group.by="Chemistry")$data[, c("features.plot", "id", "pct.exp")]
perc_exp %>%
  filter(pct.exp > MIN_PCT_GENE ) -> Acceptable_Genes

Acceptable_Genes = rownames(Acceptable_Genes)

List_total= intersect(List_Markers, Acceptable_Genes)
List_of_Lists=  lapply(List_of_Lists , function(x) intersect(x, Acceptable_Genes))





#Calculate Average expression within each cluster


PBMC$FirstClust= droplevels(PBMC$FirstClust)
PBMC= SetIdent(PBMC, value = "FirstClust")

#Try with ComplexHeatmap:
cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = List_total , slot= "data")
mat = cluster.averages$RNA
#mat.scaled= t(scale(t(mat)))
#mat.scaled = (mat-apply(mat, 1, mean) / apply(mat, 1, sd))
mat.scaled= t(apply(mat, MARGIN = 1, FUN = function(X) (X - mean(X))/sd(X)))

col_fun = circlize::colorRamp2(c(-2,-1, 0, 1,2), rev(brewer.pal(n=5,name="RdBu")))
#col_fun = circlize::colorRamp2(c(-4,-2, -1, -0.5, 0, 0.5, 1, 2, 4), rev(brewer.pal(n=9,name="RdBu")))
#col_fun = circlize::colorRamp2(c(-4,-3 , -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4), rev(brewer.pal(n=11,name="RdBu")))

SIZE_ROW = 12


p1 = Heatmap(mat.scaled[intersect(List_of_Lists[[1]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, row_names_gp = gpar(fontsize = SIZE_ROW),heatmap_legend_param = list(title="Z score", direction = "horizontal", legend_width=unit(10,"cm")) ) 
p2 = Heatmap(mat.scaled[intersect(List_of_Lists[[2]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))


#Draw the serie of heatmaps
ht_list= p1 %v% p2

draw(ht_list, heatmap_legend_side = "top")








p1 = Heatmap(mat.scaled[intersect(List_Genes1, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")) )
p2 = Heatmap(mat.scaled[intersect(List_Genes2, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")), border = TRUE )
p3 = Heatmap(mat.scaled[intersect(List_Genes3, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")), border = TRUE )
p4 = Heatmap(mat.scaled[intersect(List_Genes4, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")), border = TRUE )
p5 = Heatmap(mat.scaled[intersect(List_Genes5, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")), border = TRUE )
p6 = Heatmap(mat.scaled[intersect(List_Genes6, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")), border = TRUE )
p7 = Heatmap(mat.scaled[intersect(List_Genes7, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= rev(brewer.pal(n=4,name="RdYlBu")), border = TRUE )












cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = List_Demaria , slot= "data")



for (i in 1:6){
  cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = List_of_Lists[[i]] , slot= "data")
  print(pheatmap(cluster.averages$RNA, scale= "row", cluster_rows = TRUE, cluster_cols = FALSE, fontsize_col = 10, fontsize_row = 15))
}


cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = List_Genes1 , slot= "data")
p1 = pheatmap(cluster.averages$RNA, scale= "row", cluster_rows = TRUE, cluster_cols = FALSE, fontsize_col = 10, fontsize_row = 15)

cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = List_Genes2 , slot= "data")
p2 = pheatmap(cluster.averages$RNA, scale= "row", cluster_rows = TRUE, cluster_cols = FALSE, fontsize_col = 10, fontsize_row = 15)







pheatmap(cluster.averages$RNA, scale= "row", cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 15, fontsize_row = 10, cutree_rows = 6, gaps_row = c(10, 16, 25, 34,  39))


VlnPlot(Merged_Seurat_Rescaled, feature = "KLRC1")

head(cluster.averages[["RNA"]][, 1:5])

#Return this mean as a seuratObject

orig.levels= levels(PBMC$FirstClust)
Idents(PBMC) <- gsub(pattern = " ", replacement = "_", x = Idents(PBMC))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(PBMC) <- orig.levels
cluster.averages <- AverageExpression(PBMC, return.seurat = TRUE)
cluster.averages


DoHeatmap(cluster.averages, features = rownames(cluster.averages) , size = 3, draw.lines = FALSE, label= TRUE, group.by = "FirstClust")






DoHeatmap(PBMC, features= List_total , group.by = "FirstClust",slot= "scale.data")

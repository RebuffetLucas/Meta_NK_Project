#Create A umap OF THE USUAL SUSPECTS

library(rlist)
library(reshape2)

#Alternative with new names:
#PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")


#Load Seurat Object
PBMC_for_Metadata = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")

PBMC@meta.data = PBMC_for_Metadata@meta.data


PBMC_Norm= as.SingleCellExperiment(PBMC)
PBMC_Norm = multiBatchNorm(PBMC_Norm,  batch = PBMC_Norm$orig.ident)
PBMC_Norm= as.Seurat(PBMC_Norm)


All_genes= rownames(PBMC_Norm)
PBMC=ScaleData(PBMC_Norm, features = All_genes, do.scale = DATA_SCALE , do.center = DATA_CENTER )

PBMC$FirstClust = factor(PBMC$FirstClust, levels= ORDER_CLUST_LEGEND)




FeaturePlot(PBMC_for_Metadata, features= "ZNF683", pt.size  =0.4, max.cutoff = 'q50')  + scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu")))
VlnPlot(PBMC, features= "TRIB1", group.by = "SecondClust") 

VlnPlot(Seurat_NK2 , features= "TRIB1") 

VlnPlot(PBMC, features= "NCR1") 


#Get the list of genes present in at list X percent of the cells


MIN_PCT_GENE= 5


KIR_inthenames = rownames(PBMC)[grepl("KIR", rownames(PBMC))]


#NK2

PBMC= SetIdent(PBMC, value = "FirstClust")

PBMC$FirstClust= droplevels(PBMC$FirstClust)



List_BrightvsDim = c("XCL1", "XCL2", "GZMK" , "FGFBP2", "GZMB", "KLRC1", "FCGR3A")



List_Genes1 = c( "CXCL10", "CXCL13", "CSF1", "XCL1", "XCL2", "TNF", "IFNG", "CCL3", "CCL4", "CCL5", "CSF2", "FLT3LG", "IL1B" , "IL16", "IL17C","IL18","IL32", "CCL4L2", "CCL3L3" )

List_Genes2 = c("CXCR6", "CXCR2", "CXCR3", "CXCR1", "S1PR4", "S1PR1", "CXCR4","S1PR5", "CCR1", "CCR4","CCR5", "CCR7", "CX3CR1")

#Activating R

List_Genes3 = c("SLAMF6", "KLRC2","SLAMF7","CD160", "NCR1", "NCR2", "NCR3",  "KLRK1")

#Inhibiting R
List_Genes4 = c("HAVCR2",  "LAG3", "TIGIT", "PDCD1", "KLRC1", "KIR3DL3" ,"KIR2DL3", "KIR2DL1" ,"KIR2DL4", "KIR3DL1", "KIR3DL2", "CD300A", "SIGLEC7","SIGLEC9", "KLRB1")


List_Genes5 = c("IL12RB1", "IL2RB1", "IL2RG", "IL21R", "IL18R1", "IL4R",  "IL12RB2", "IL18RAP", "IL2RA","IL2RB", "IL15RA", "IL10RA", "IL10RB", "IL7R", "IL21R", "IL18R1", "TGFBR1", "TGFBR3", "TGFBR2")


List_Genes6 = c("TNFSF10", "FASLG", "PRF1", "GZMB" , "NKG7", "GZMH", "GZMA", "GSDMB" ,  "GSDMD", "GZMK" )

List_Genes7 =c("TNFRSF9", "CD69", "TNFRSF18", "TNFSF4", "RGS1")

List_Genes8 = c("NCAM1","ITGAX", "ITGAM", "B3GAT1", "FCGR3A")


List_Markers = list.append(List_Genes1, List_Genes2, List_Genes3, List_Genes4, List_Genes5, List_Genes6, List_Genes7, List_Genes8)
List_of_Lists= list(List_Genes1, List_Genes2, List_Genes3, List_Genes4, List_Genes5, List_Genes6 , List_Genes7, List_Genes8)

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
p3 = Heatmap(mat.scaled[intersect(List_of_Lists[[3]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))
p4 = Heatmap(mat.scaled[intersect(List_of_Lists[[4]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))
p5 = Heatmap(mat.scaled[intersect(List_of_Lists[[5]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))
p6 = Heatmap(mat.scaled[intersect(List_of_Lists[[6]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))
p7 = Heatmap(mat.scaled[intersect(List_of_Lists[[7]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))
p8 = Heatmap(mat.scaled[intersect(List_of_Lists[[8]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))


#Draw the serie of heatmaps
ht_list= p1 %v% p2 %v% p3 %v% p4 %v% p5 %v% p6 %v% p7  %v% p8

draw(ht_list, heatmap_legend_side = "top")
p14 = draw(ht_list, heatmap_legend_side = "top")

p14

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/HeatMap_Usualsupects_V2.png", width = 18, height = 30,  units = "cm", res=600 )
p14
dev.off()

#Fig3a
pdf(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig3a.pdf", width = 8, height = 15)
p14
dev.off()











VlnPlot(PBMC, features= c("GZMB"))


DimPlot(PBMC_for_Metadata, group.by = "nkg2c")



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

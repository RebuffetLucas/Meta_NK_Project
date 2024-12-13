#Draw a heatmap of common NK surface prots (list from Marie_Laure)

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


List_Genes1 = c(  "NCAM1", "CD7", "CD2", "FCGR3A", "KIR2DL1", "KIR2DL2", "KIR2DS1", "KIR2DS2", "KIR3DL1", "KIR3DL2","KLRC1","KLRC2", "KLRK1", "CD226", "CD244", "NCR3", "KLRF1", "NCR2" , "NCR1" )
List_Genes2 = c("B3GAT1", "PTPRC", "FCER1G", "ITGA1", "TNFRSF9", "LAMP1", "ITGAE","CD69", "IL2RA", "IL7R", "IL2RB", "IL2RG", "IL18R", "LILRB2", "LILRB1", "SIGLEC7", "CXCR3", "CCR5")

#intersect(All_genes , c("B3GAT1", "CD244", "CD25", "IL2RB", "IL2RG", "SIGLEC7"))

List_Markers = list.append(List_Genes1, List_Genes2)
List_of_Lists= list(List_Genes1, List_Genes2)

List_total= intersect(List_Markers, rownames(PBMC))


#Keep Genes expressed in more than X percent cells
perc_exp <- DotPlot(PBMC, features=List_total, group.by="FirstClust")$data[, c("features.plot", "id", "pct.exp")]

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







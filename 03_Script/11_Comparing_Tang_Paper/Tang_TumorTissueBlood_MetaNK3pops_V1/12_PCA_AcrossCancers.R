library(ade4)
library(magrittr)
#install.packages("factoextra") 
library(factoextra)

#Load data
data.ref_metaNK = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_TumorTissueBlood_MetaNK3pops_V1/Ref_Integrated.rds")
data.query_blood  = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_Blood_MetaNK3pops_V1/data.query_all_cells_BloodTang.rds")
data.query_tumor  = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_TumorTissueBlood_MetaNK3pops_V1/data.query_all_cells_Tumor_Normal.rds")

#Keep only tumor, not peripheral NK
data.query_tumor= subset(data.query_tumor, subset= meta_tissue == "Tumor")

#Have a look at the data + add a triple identity slot
data.ref_metaNK$meta_tissue = "Blood"
table(data.ref_metaNK$meta_histology)
table(data.ref_metaNK$FirstClust)

table(data.query_blood$meta_tissue)
table(data.query_blood$meta_histology)
data.query_blood$FirstClust = data.query_blood$predicted.id
table(data.query_blood$FirstClust)


table(data.query_tumor$meta_tissue)
table(data.query_tumor$meta_histology)
data.query_tumor$FirstClust = data.query_tumor$predicted.id
table(data.query_tumor$FirstClust)


#Merge them together
data_merged = merge(x=data.ref_metaNK , y =c(data.query_blood , data.query_tumor) , add.cell.ids = c("ref", "Tang_Blood", "Tang_Tumor") )
DefaultAssay(data_merged) = "RNA"

#Abreviation of cancer
data_merged$meta_histology_relabel  = as.factor(data_merged$meta_histology)
levels(data_merged$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(data_merged$meta_histology_relabel ))

#double ID
data_merged$tissue_Histo = paste0(data_merged$meta_tissue , data_merged$meta_histology_relabel)

#Triple ID
data_merged$Tissue_Histo_FirstClust = paste0(data_merged$meta_tissue , "_", data_merged$meta_histology_relabel,"_", data_merged$FirstClust)


#Find Variable features

List_Tissue_Histo = SplitObject(data_merged, split.by = "tissue_Histo")

for (i in 1:length(List_Tissue_Histo)) {
  List_Tissue_Histo[[i]] <- NormalizeData(List_Tissue_Histo[[i]], verbose = FALSE) #Maybe we need to change to MultiBatchNorm Here
  List_Tissue_Histo[[i]] <- FindVariableFeatures(List_Tissue_Histo[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


#Variable Features
FEATURES_RETAINED = SelectIntegrationFeatures(
  List_Tissue_Histo,
  nfeatures = 2000,
  assay = rep("RNA", length(List_Tissue_Histo)) ,
  verbose = TRUE
)



#Merged_Seurat_Rescaled = merge(List_Tissue_Histo[[1]] , y = List_Tissue_Histo[-1])
Merged_Seurat_Rescaled = merge(List_Tissue_Histo[[2]] , y = c( List_Tissue_Histo[[3]], List_Tissue_Histo[[4]], List_Tissue_Histo[[5]])) #Just a test here
Merged_Seurat_Rescaled = SetIdent(Merged_Seurat_Rescaled, value= "tissue_Histo")


if (DO_SUBSET_FOR_PCA == TRUE){
Merged_Seurat_Rescaled = subset(Merged_Seurat_Rescaled, downsample = NUMBER_CELLS_PERCATEGORY_PCA )
}

if (DO_SUBSET_TISSUE_FOR_PCA == TRUE){
  Merged_Seurat_Rescaled = subset(Merged_Seurat_Rescaled, subset = "meta_tissue"== TISSUE_TO_SUBSET_FOR_PCA )
  Merged_Seurat_Rescaled=ScaleData(Merged_Seurat_Rescaled, features =rownames(Merged_Seurat_Rescaled) , do.scale = DATA_SCALE , do.center = DATA_CENTER )
}

Merged_Seurat_Rescaled=ScaleData(Merged_Seurat_Rescaled, features =rownames(Merged_Seurat_Rescaled) , do.scale = DATA_SCALE , do.center = DATA_CENTER )


#Look at the mean
daframe_Average = AverageExpression(
  Merged_Seurat_Rescaled,
  assays = "RNA",
  features = FEATURES_RETAINED,
  return.seurat = FALSE,
  group.by = "Tissue_Histo_FirstClust",
  add.ident = NULL,
  layer = "data",
  slot= "scale.data",
  verbose = TRUE
)

mat = daframe_Average$RNA

mat= t(mat)

#Run PCA with ade4

res.pca <- dudi.pca(mat,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5 ,          # Nombre d'axes gardés
                    center= TRUE ,
                    scale= TRUE
)


#Basic
fviz_eig(res.pca)


mat[1:5,1:5]

df= as.data.frame( mat )
df$group = str_sub(rownames(df),-3,-1)
df$group = as.factor(df$group)
  

#PCA viz data colr
p2 = fviz_pca_ind(res.pca,
             col.ind = df$group, 
             palette = palette3,
             repel = TRUE ,
             labelsize= 7,
             pointsize= 2,
             mean.point= FALSE
)  + theme(text = element_text(size = 20))  +   labs(title ="PCA", x = "PC1", y = "PC2")       



png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_All_Cancers", "/PCA_No_Elypse.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p2)
dev.off()




#With Elypses

p1 = fviz_pca_ind(res.pca,
             #col.ind = "coord", 
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             habillage=df$group,
             axes= c(1,2),
             cols.ind = palette3,
             invisible='quali',
             repel = TRUE ,
             labelsize= 8,
             pointsize= 3,
             addEllipses= TRUE,
             #ellipse.level=1.8,
             legend= "none"
)  + theme(text = element_text(size = 20)) #+ ggforce::geom_mark_ellipse(aes(fill = df$group ,  color = df$group)) +  coord_equal()

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_All_Cancers", "/PCA_With_Elypse.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()












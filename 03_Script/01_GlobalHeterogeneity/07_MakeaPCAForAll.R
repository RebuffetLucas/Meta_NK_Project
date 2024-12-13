library(ade4)
library(magrittr)
library(factoextra)

install.packages("ggforce")
remotes::install_version("ggforce" , "0.4.0")
#devtools::install_github("thomasp85/ggforce")

#Make a PCA for all clusters
PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

#Keeping only the genes conserved across samples
minrow= rownames(PBMC)
PBMC = subset(PBMC , features = minrow)
PBMC= SetIdent(PBMC, value = "FirstClust")

#ONLY VariABLE FEATURES:
VF1= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VariableFeatures.rds")

FEATURES_TO_USE= VF1

#Normalization

sce_Merged= as.SingleCellExperiment(PBMC)
sce_Merged_Rescaled = multiBatchNorm(sce_Merged,  batch = sce_Merged$orig.ident)
Merged_Seurat_Rescaled= as.Seurat(sce_Merged_Rescaled)



PBMC=ScaleData(Merged_Seurat_Rescaled, features = minrow, do.scale = DATA_SCALE , do.center = DATA_CENTER )





#Try with ComplexHeatmap:
cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = FEATURES_TO_USE , slot= "scale.data")
mat = cluster.averages$RNA

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

df= as.data.frame( mat[,1] )
df$group=c("NK1","NK3","NK1" ,"NK1", "NK2" , "NK2")
df$group = as.factor(df$group)

fviz_pca_ind(res.pca,
             col.ind = "coord", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ,
             labelsize= 7,
             pointsize= 2
)  + theme(text = element_text(size = 20))  +   labs(title ="PCA", x = "PC1", y = "PC2")       




#Colored by main pops HAS TO BE FIXED

fviz_pca_ind(res.pca,
             col.ind = "coord", 
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ,
             labelsize= 8,
             habillage=df$group,
             pointsize= 3,
             label="none"
)  + theme(text = element_text(size = 20))             






#Try to put ellipses HAS TO BE FIXED
df$group  = factor(df$group, levels = c("NK1", "NK3", "NK2"))
fviz_pca_ind(res.pca,
             #col.ind = "coord", 
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             habillage=df$group,
             axes= c(1,2),
             invisible='quali',
             repel = TRUE ,
             labelsize= 8,
             pointsize= 3,
             addEllipses= TRUE,
             #ellipse.level=1.8,
             legend= "none"
)  + theme(text = element_text(size = 20)) + ggforce::geom_mark_ellipse(aes(fill = df$group ,  color = df$group)) +  coord_equal()



df$group  = factor(df$group, levels = c("NK1", "NK3", "NK2"))
p_PCA = fviz_pca_ind(res.pca,
             #col.ind = "coord", 
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             habillage=df$group,
             axes= c(1,2),
             invisible='quali',
             repel = TRUE ,
             labelsize= 8,
             pointsize= 3,
             addEllipses= TRUE,
             #ellipse.level=1.8,
             legend= "none"
)  + theme(text = element_text(size = 20))  +  coord_equal()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/PCA_Without_Elipse.png", width = 20, height = 20,  units = "cm", res=600 )
p_PCA
dev.off()




#Try with ComplexHeatmap:
cluster.averages <- AverageExpression(PBMC, group.by = c("SecondClust","orig.ident"), features = minrow , slot= "data")
mat = cluster.averages$RNA
mat= t(mat)










#Optimization of the granularity of the first clustering
#Using Sab's method to assess the purity of the clusters

install.packages("clustree")

devtools::install_github("PaulingLiu/ROGUE")

library(clustree)
library(ggraph)

library(ROGUE)
library(ggplot2)
library(tidyverse)

#PBMC_Sab= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK4/05_Output/01_GlobalHeterogeneity/PBMC_All.rds")

PBMC_Sab = Merged_Seurat_Rescaled
#Check the object
DimPlot(PBMC_Sab, reduction = "UMAP")
FeaturePlot(PBMC_Sab, reduction = "UMAP", feature= c("XCL1","GZMB", "GZMK", "CXCR4"))


# Remove previous clustering results
PBMC_Sab$seurat_clusters <- NULL
PBMC_Sab@meta.data <- PBMC_Sab@meta.data %>% dplyr::select(-matches("^RNA_snn_res."))

#DefaultAssay(PBMC_Sab ) <- "integrated"


# The first step aims to compute cells clusters using a graph-based clustering approach

PBMC_Sab <- FindNeighbors(object = PBMC_Sab, 
                              reduction = "harmony",
                              dims = 1:30, 
                              verbose = FALSE,
                              force.recalc = TRUE )


for(i in seq(0.5, 1.4 , 0.1)){
  PBMC_Sab <- FindClusters(object = PBMC_Sab, 
                               resolution = i, 
                               verbose = FALSE,
                               random.seed = SEED)
}




# Choose cluster optimal resolution : visualize cluster stability per res--------------------------

# Compute clustree
clustree.out <- clustree(PBMC_Sab ,prefix = "RNA_snn_res.", layout="sugiyama", node_colour = "sc3_stability", node_text_colour = "white")
print(clustree.out)

# Compute mean stability per res
stability.df <- clustree.out[["data"]] %>% dplyr::group_by(RNA_snn_res.) %>% dplyr::summarise(mean.stability=mean(sc3_stability))


# Compute difference between 2 successive point of resolution
stability.df$difference <- NA
for(i in 1:(nrow(stability.df)-1)){
  stability.df$difference[i] <- stability.df$mean.stability[i+1] - stability.df$mean.stability[i]
}

# Store the optimal res. in seurat meta data
opt.res <- stability.df %>% filter(difference < 0 | is.na(difference)) %>% filter(mean.stability==max(mean.stability)) %>% mutate(RNA_snn_res.= paste0("RNA_snn_res.",RNA_snn_res.)) %>% pull(RNA_snn_res.)


PBMC_Sab[["integrated.optimal.resolution"]] <- PBMC_Sab[[opt.res[1]]]
PBMC_Sab[["integrated.optimal.resolution"]] <- paste0("integrated.Cl.",PBMC_Sab$integrated.optimal.resolution)



# Visualize stability per res
stability.df$RNA_snn_res. <- paste0("RNA_snn_res.",stability.df$RNA_snn_res.)
ggplot(stability.df, aes(x=RNA_snn_res. , y=mean.stability, group = 1)) + 
  geom_point(color = "steelblue") + 
  geom_line(color = "steelblue") + 
  labs(x="")+
  ylim(0,(max(stability.df$mean.stability)+0.001)) + 
  geom_vline(xintercept = opt.res[1], linetype="dashed", color = "orange" ) + 
  ggtitle("Optimal resolution", subtitle = "Mean sc3 Stability Method") + 
  theme_classic() + 
  theme( plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size = 14), 
         plot.subtitle = element_text(lineheight=.8, face="bold", hjust = 0.5, size = 12),
         axis.title.y = element_text(lineheight=.8, hjust = 0.5, size = 12),
         axis.text = element_text(lineheight=.8, size = 12, angle = 45,  vjust = 1, hjust = 1))


DimPlot(PBMC_Sab, reduction = "UMAP", group.by = "RNA_snn_res.0.7", label = TRUE, label.size = )

DimPlot(PBMC_Sab, reduction = "UMAP", group.by = "RNA_snn_res.0.6")

DimPlot(PBMC_Sab, reduction = "umap", group.by = "RNA_snn_res.1.2", pt.size = 1.8 )

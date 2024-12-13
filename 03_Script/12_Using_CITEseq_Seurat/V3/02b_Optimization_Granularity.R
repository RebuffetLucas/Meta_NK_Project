#Optimization of the granularity of the first clustering
#Using Sabrina's method to assess the purity of the clusters

install.packages("clustree")
install.packages("tidyverse")

devtools::install_github("PaulingLiu/ROGUE")

library(clustree)
library(ggraph)

library(ROGUE)
library(ggplot2)
library(tidyverse)

#PBMC_Sabrina= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK4/05_Output/01_GlobalHeterogeneity/PBMC_All.rds")

PBMC_Sabrina = Seurat_NK_Final
#Check the object
DimPlot(PBMC_Sabrina, reduction = "wnn.umap")
FeaturePlot(PBMC_Sabrina, reduction = "wnn.umap", feature= c("XCL1","GZMB", "GZMK", "CXCR4"))


# Remove previous clustering results
PBMC_Sabrina$seurat_clusters <- NULL
PBMC_Sabrina@meta.data <- PBMC_Sabrina@meta.data %>% dplyr::select(-matches("^wsnn_res."))

#DefaultAssay(PBMC_Sabrina ) <- "integrated"


# The first step aims to compute cells clusters using a graph-based clustering approach

#PBMC_Sabrina <- FindNeighbors(object = PBMC_Sabrina, reduction = "harmony", dims = 1:30,  verbose = FALSE, force.recalc = TRUE )


for(i in seq(0.2, 0.8 , 0.05)){
  PBMC_Sabrina <- FindClusters(object = PBMC_Sabrina,
                               graph.name = "wsnn",
                               resolution = i,
                               algorithm = 3,
                               verbose = TRUE,
                               random.seed = SEED)
}




# Choose cluster optimal resolution : visualize cluster stability per res--------------------------

# Compute clustree
clustree.out <- clustree(PBMC_Sabrina ,prefix = "wsnn_res.", layout="sugiyama", node_colour = "sc3_stability", node_text_colour = "white")
print(clustree.out)

# Compute mean stability per res
stability.df <- clustree.out[["data"]] %>% dplyr::group_by(wsnn_res.) %>% dplyr::summarise(mean.stability=mean(sc3_stability))


# Compute difference between 2 successive point of resolution
stability.df$difference <- NA
for(i in 1:(nrow(stability.df)-1)){
  stability.df$difference[i] <- stability.df$mean.stability[i+1] - stability.df$mean.stability[i]
}

# Store the optimal res. in seurat meta data
opt.res <- stability.df %>% filter(difference < 0 | is.na(difference)) %>% filter(mean.stability==max(mean.stability)) %>% mutate(wsnn_res.= paste0("wsnn_res.",wsnn_res.)) %>% pull(wsnn_res.)


PBMC_Sabrina[["integrated.optimal.resolution"]] <- PBMC_Sabrina[[opt.res[1]]]
PBMC_Sabrina[["integrated.optimal.resolution"]] <- paste0("integrated.Cl.",PBMC_Sabrina$integrated.optimal.resolution)



# Visualize stability per res
stability.df$wsnn_res. <- paste0("wsnn_res.",stability.df$wsnn_res.)
ggplot(stability.df, aes(x=wsnn_res. , y=mean.stability, group = 1)) + 
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


DimPlot(PBMC_Sabrina, reduction = "wnn.umap", group.by = "wsnn_res.0.1", label = TRUE, label.size = )

DimPlot(PBMC_Sabrina, reduction = "wnn.umap", group.by = "wsnn_res.0.6")

DimPlot(PBMC_Sabrina, reduction = "wnn.umap", group.by = "wsnn_res.1.2", pt.size = 1.8 )

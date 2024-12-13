#This script run the trajectory analysis for the entire subsets of Dim together

#For Dim

# Perform the transcriptionnal trajectories analysis (Monocle 3) on NK Cells
library(monocle3)
library(Seurat)
library(SeuratData)
library(viridis)

# Version of seurat Wrapper that actually works
# remotes::install_github('satijalab/seurat-wrappers', ref = "b8feea013e7e19a46e935684b510399ffe0b6740" #
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

library(cerebroApp)



#Start from suerat object

PBMC= NK_Seurat

NK_Seurat= PBMC

DimPlot(NK_Seurat)

NK_Seurat  =subset(PBMC, idents= CLUSTERS_TO_REMOVE, invert= TRUE)
cds= as.cell_data_set(NK_Seurat)
cds= estimate_size_factors(cds)


#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Assign cluster informations
list.cluster <- NK_Seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster


#Plot before trajectory
partition.before.traj <-plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
partition.before.traj


cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


#Learn Trajectory
#cds <- learn_graph(cds, use_partition = F, learn_graph_control  = list(ncenter=280), close_loop = FALSE) #Best so far

#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE, learn_graph_control  = list(minimal_branch_len=25))

cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE)


cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE,  learn_graph_control=LIST_PARAM_CONTROL)



#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE)


#cds <- learn_graph(cds,  use_partition = FALSE,learn_graph_control = list( nn.k=30), close_loop = FALSE)

#cds <- learn_graph(cds, learn_graph_control = list( nn.k=20), use_partition = FALSE, close_loop = FALSE)

#cds <- learn_graph(cds, learn_graph_control  = list( nn.k=15), close_loop = FALSE)
#cds <- learn_graph(cds, use_partition = F, close_loop = FALSE)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 2)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)

p1 = plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)

#Keep on with the Dynamic UMAP

#Order cells in PseudoTime
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

p2 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)

p3 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE , label_cell_groups= TRUE,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5, 
           group_label_size =  15, label_principal_points = F, trajectory_graph_segment_size = 1.8, trajectory_graph_color= "blue")

p4 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE , label_cell_groups= TRUE,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5, 
                group_label_size =  15, label_principal_points = F, trajectory_graph_segment_size = 1.8, trajectory_graph_color= "red")



#Filter the most significant
   #Best of the best
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

  #Top50
modulated_genes %>% dplyr::arrange(q_value, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=50) -> top50

row.names(top50)

  #Top100
modulated_genes %>% dplyr::arrange(q_value<0.05, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(!grepl("RPS", rownames(.))) %>%
  dplyr::filter(!grepl("RPL", rownames(.))) %>%
  dplyr::filter(!grepl("MT-", rownames(.))) %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=150) -> top100

row.names(top100)


genes=row.names(top100)


#Order them by pseudotime
pt.pseudotime = as.data.frame(pseudotime(cds)[order(pseudotime(cds))])
pt.pseudotime = t(pt.pseudotime)

#pt.pseudotime= as.matrix(pt.pseudotime)
#pt.pseudotime = rbind(rep(0,  length(pt.pseudotime)), pt.pseudotime)


pt.pseudotime = t(as.matrix(pseudotime(cds)[order(pseudotime(cds))]))


htpseudo <- Heatmap(
  pt.pseudotime,
  name                         = "pseudotime",
  col                          = plasma(11, alpha=1, begin = 0, end=1, direction = 1) ,
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  km = 1,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =unit(0.5, "cm") )


print(htpseudo)

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]


#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

HEIGHT = unit(15, "cm")
WIDTH = unit(12, "cm")
POLICE= gpar(fontsize = 3)



NUMBER_CLUSTS = 10


#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  km = NUMBER_CLUSTS,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = TRUE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH,
  row_title = NULL)

ht_list = htpseudo %v% htkm
draw(ht_list)


FILE = paste0("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/Dynamic_Hetamap_V2/DynamicHeatmapTop150_NKintand1", NUMBER_CLUSTS, "clusters",".png")

png(file= FILE, width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list)
dev.off()



#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = TRUE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH)



ht_list2 = htpseudo %v% hthc
draw(ht_list2)

FILE = paste0("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/Dynamic_Hetamap_V2/DynamicHeatmapTop150_NKintand1", "WardMETHOD",".png")

png(file= FILE, width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list2)
dev.off()


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajOntheUMAP.png", width = 30, height =25, units= "cm" ,  res=600 )
p1
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajPseudotime.png", width = 30, height =25, units= "cm" ,  res=600 )
p2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajPseudotimeblue.png", width = 30, height =25, units= "cm" ,  res=600 )
p3
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajPseudotimered.png", width = 30, height =25, units= "cm" ,  res=600 )
p4
dev.off()





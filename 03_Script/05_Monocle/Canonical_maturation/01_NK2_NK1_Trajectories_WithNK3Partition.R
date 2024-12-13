#This script run the trajectory analysis for the entire subsets of Dim together

#For Dim

# Perform the transcriptionnal trajectories analysis (Monocle 3) on NK Cells
library(monocle3)
library(Seurat)
library(SeuratData)

# Version of seurat Wrapper that actually works
# remotes::install_github('satijalab/seurat-wrappers', ref = "b8feea013e7e19a46e935684b510399ffe0b6740" #
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

#Start from suerat object

PBMC= NK_Seurat

NK_Seurat= PBMC


cds= as.cell_data_set(NK_Seurat)
cds= estimate_size_factors(cds)


#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions


#Assign cluster informations
list.cluster <- NK_Seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

#Assign partition  2 to NK3

NK3_cells =  names(list.cluster[which(list.cluster==c("NK3"))])




levels(recreate.partitions) = c(1,2)
for(name in NK3_cells) {
  recreate.partitions[[name]] <- factor("2", levels = c("1", "2"))
}

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions


#Plot before trajectory
partition.before.traj <-plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
partition.before.traj


cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


#Learn Trajectory
#cds <- learn_graph(cds, use_partition = F, learn_graph_control  = list(ncenter=280), close_loop = FALSE) #Best so far

cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE, learn_graph_control  = list(minimal_branch_len=30))
cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE, learn_graph_control=list(ncenter=400))

#cds <- learn_graph(cds,  use_partition = FALSE,learn_graph_control = list( nn.k=30), close_loop = FALSE)

#cds <- learn_graph(cds, learn_graph_control = list( nn.k=20), use_partition = FALSE, close_loop = FALSE)

#cds <- learn_graph(cds, learn_graph_control  = list( nn.k=15), close_loop = FALSE)
#cds <- learn_graph(cds, use_partition = F, close_loop = FALSE)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 2)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)



#Stop here to pursue with Branch analysis




#Order cells in PseudoTime
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

#Visualize Pseutotime
ggplot(data.pseudo, aes(monocle3_pseudotime, Secondclust, fill = Secondclust )) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Secondclust, monocle3_pseudotime), fill = Secondclust)) + geom_boxplot()

#Find Genes that change as a function of pseudotime

#deg <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4) #Long step
deg= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/03_Transcriptionnal_Trajectories/ForDim/deg_NK1.rds")
deg %>% dplyr::arrange(q_value) %>%
  dplyr::filter(status == "OK") %>%
  head(n=10) -> top10

#Add pseudotime values into the seuratobject
NK_Seurat$pseudotime <- pseudotime(cds)
FeaturePlot(NK_Seurat, features = "pseudotime")

#saveRDS(NK_Seurat, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/03_Transcriptionnal_Trajectories/ForDim/NK1withPseudotime.rds")

#Have a look at genes of interest
RidgePlot(NK_Seurat, features = rownames(top10), sort = T, idents = c("NK_1A", "NK_1D", "NK_1B"))

#genes of interest
fData(cds)$gene_short_name <- rownames(fData(cds))
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("PRDM2", "RUNX3", "TBX21"))) 
cds_subset <- cds[my_genes,]

plot_genes_in_pseudotime(cds_subset, color_cells_by = "Secondclust", min_expr = 0.4 )

saveRDS(cds, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/05_Output/03_Transcriptionnal_Trajectories/ForDim/NK1_cds_Object.rds" )




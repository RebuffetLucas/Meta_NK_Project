#Try to extract trajectory and inject into complete object


# WIth the function that does not wor, only compatible with old monocle
#remotes::install_github("romanhaa/cerebroApp")

#NK_Seurat2 = extractMonocleTrajectory(cds, NK_Seurat , 
#                                      
#                                      column_state = 'State',
#                                      column_pseudotime = 'Pseudotime' )

#Extract pseudotime
traj.coord<- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

#Extract traj points:
traj.plot <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
                        label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                        group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)




point.data <- ggplot_build(traj.plot)[["plot"]][["data"]]

p = ggplot_build(traj.plot)[["data"]][[4]]



#Extract coordinates for the branches
branches.data = ggplot_build(traj.plot)[["data"]][[3]]


#Build new cds2 with all the cells
cds2= as.cell_data_set(PBMC)
cds2= estimate_size_factors(cds2)


#Assign partitions
recreate.partitions2 <- c(rep(1, length(cds2@colData@rownames)))
names(recreate.partitions2) <- cds2@colData@rownames
recreate.partitions2 <- as.factor(recreate.partitions2)
recreate.partitions2


#Assign cluster informations
list.cluster <- PBMC@active.ident
cds2@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

#Assign partition  2 to NK3

#NK3_cells =  names(list.cluster[which(list.cluster=="NK3")])

NK3_cells =  names(list.cluster[which(list.cluster%in% CLUSTERS_TO_REMOVE) ])



levels(recreate.partitions2) = c(1,2)
for(name in NK3_cells) {
  recreate.partitions2[[name]] <- factor("2", levels = c("1", "2"))
}

cds2@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions2

#Plot before trajectory
partition.before.traj <-plot_cells(cds2, color_cells_by = "partition", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
partition.before.traj


cluster.before.traj <-plot_cells(cds2, color_cells_by = "cluster", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj

#Before injecting
plot_cells(cds2, color_cells_by = "cluster", label_groups_by_cluster = F , group_label_size= 0,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 0, label_principal_points = T, trajectory_graph_segment_size = 2)

#Learn Trajectory

cds2= learn_graph(cds2)

cds2 = order_cells(cds2, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

#Order cells

plot_cells(cds2, color_cells_by = "cluster", label_groups_by_cluster = F , group_label_size= 0,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 0, label_principal_points = T, trajectory_graph_segment_size = 2)


plot_cells(cds2, color_cells_by = "pseudotime", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)


#have a look
cds3 = cds2


cds3@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][names(traj.coord)] = traj.coord

p1_v2 = plot_cells(cds3, color_cells_by = "cluster", label_groups_by_cluster = FALSE , label_cell_groups = FALSE, 
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2, show_trajectory_graph = FALSE) +  geom_segment(data = branches.data,aes(x=x, y=y, xend=xend, yend=yend), color="grey28", size=2, linetype="solid", alpha=0.8) +  labs(title="Branches Plot") +  theme_classic() +   scale_color_manual(values=palette)
p1_v2 

p2_v2 = plot_cells(cds3, color_cells_by = "pseudotime", label_groups_by_cluster = F , label_cell_groups = FALSE,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2, show_trajectory_graph = FALSE) +  geom_segment(data = branches.data,aes(x=x, y=y, xend=xend, yend=yend), color="grey28", size=2, linetype="solid", alpha=0.8) +  labs(title="Branches Plot") +  theme_classic()

p2_v2


p3_v2 = plot_cells(cds3, color_cells_by = "pseudotime", label_groups_by_cluster = F , label_cell_groups = FALSE,
                   label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                   graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2, show_trajectory_graph = FALSE) +  geom_segment(data = branches.data,aes(x=x, y=y, xend=xend, yend=yend), color="blue", size=2, linetype="solid", alpha=0.8) +  labs(title="Branches Plot") +  theme_classic()

p3_v2


p4_v2 = plot_cells(cds3, color_cells_by = "pseudotime", label_groups_by_cluster = F , label_cell_groups = FALSE,
                   label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                   graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2, show_trajectory_graph = FALSE) +  geom_segment(data = branches.data,aes(x=x, y=y, xend=xend, yend=yend), color="red", size=2, linetype="solid", alpha=0.8) +  labs(title="Branches Plot") +  theme_classic()

p4_v2




#Save in HD

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajOntheUMAP_V2_From_NKint.png", width = 30, height =25, units= "cm" ,  res=600 )
p1_v2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajPseudotime_V2_From_NKint.png", width = 30, height =25, units= "cm" ,  res=600 )
p2_v2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajPseudotimeblue_V2.png", width = 30, height =25, units= "cm" ,  res=600 )
p3_v2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/05_Monocle/TrajPseudotimered_2.png", width = 30, height =25, units= "cm" ,  res=600 )
p4_v2
dev.off()



#Fig4f
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig4f.pdf",
  plot = p1_v2,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 30,
  height = 25,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#Fig4g
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig4g.pdf",
  plot = p2_v2,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 30,
  height = 25,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)









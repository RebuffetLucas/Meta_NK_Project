#Optimize first trajectory research


#Optimize both ncenter and K:
ncenter = seq(250,300,10)
k = seq(15,50,5)

for (i in ncenter){
  for (j in k){
    
    print(paste0("ncenter", i, "K", j))
    
    cds <- learn_graph(cds, learn_graph_control  = list(ncenter= i , nn.k = j), close_loop = FALSE)
    
    
    plot(plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
               label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
               graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 1.25)    + ggtitle(paste0("ncenter", i, "K", j)))
     
  }
}


#Optimize only K:
for (j in k){
  
  print(paste0("ncenter", "NONE", "K", j))
  
  cds <- learn_graph(cds, learn_graph_control  = list( nn.k = j), close_loop = FALSE)
  
  plot(plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
                  label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                  graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 1.25)    + ggtitle(paste0("ncenter", "NONE", "K", j)))
  
}

## @knitr DiagnosticAfterLabelTransfer

cat("# Diagnostic After Label Transfer {.tabset .tabset-fade} \n\n")

cat("## Meta NK data raw \n\n")


  #Global UMAP
p01 = DimPlot(PBMC_Meta_NK_V2, cols = palette2) 


cat("\n\n")
print(p01)
cat("\n\n")




cat("## Vizualization of the reference {.tabset .tabset-fade} \n\n")

  #Have a look at the reference built with the 4 Dataset of Meta-NK

p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "Dataset")
p2 <- DimPlot(data_integrated, reduction = "umap", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette2) + NoLegend()

cat("\n\n")
print( p1 + p2 )

cat("\n\n")

p4 <- DimPlot(data_integrated, reduction = "umap", split.by = "Dataset", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette2) + NoLegend()

print(p4)
cat("\n\n")


cat("## Assessing quality of Label Transfer {.tabset .tabset-fade} \n\n")

cat("### Look at the common markers for the different populations {.tabset .tabset-fade} \n\n")

data.query = SetIdent(data.query, value = "predicted.id")


cat("#### Main pops based on V2 \n\n")

  #Main pops
Markers_NK_123 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds")
  
  #ViolinPlot
data.query_loop = data.query
  
for (i in names(Markers_NK_123)){
  data.query_loop = AddModuleScore(data.query_loop, features = as.list(Markers_NK_123[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln1 = VlnPlot(data.query_loop, features = paste0(names(Markers_NK_123),"1") , group.by = "predicted.id", cols= palette2, pt.size = 0)   

cat("\n\n")
print(p_Vln1)

cat("\n\n")


  #Direct Dot Plot
p5 = DotPlot(data.query_loop, features = Markers_NK_123$NK1, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) 
p6 = DotPlot(data.query_loop, features = Markers_NK_123$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p7 = DotPlot(data.query_loop, features = Markers_NK_123$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p5)
cat("\n\n")
print(p6)
cat("\n\n")
print(p7)
cat("\n\n")

cat("#### Main pops based on CITEseq \n\n")

#Main pops
Markers_NK_123_CITEseq = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds")

Markers_NK_123_CITEseq %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All

top_All %>% filter(cluster== "NK_1") -> Top_NK1
top_All %>% filter(cluster== "NK_2") -> Top_NK2
top_All %>% filter(cluster== "NK_3") -> Top_NK3


#ViolinPlot
data.query_loop = data.query
data.query_loop = AddModuleScore(data.query_loop, features = list(Top_NK1$gene), pool= NULL ,name= "NK1" , seed=19)
data.query_loop = AddModuleScore(data.query_loop, features = list(Top_NK2$gene), pool= NULL ,name= "NK2" , seed=19)
data.query_loop = AddModuleScore(data.query_loop, features = list(Top_NK3$gene), pool= NULL ,name= "NK3" , seed=19)


p_Vln2 = VlnPlot(data.query_loop, features = c("NK11", "NK21", "NK31") , group.by = "predicted.id", cols= palette2, pt.size = 0)   

cat("\n\n")
print(p_Vln2)
cat("\n\n")



p8 = DotPlot(data.query, features = Top_NK1$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) 
p9 = DotPlot(data.query, features = Top_NK2$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p10 = DotPlot(data.query, features = Top_NK3$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p8)
cat("\n\n")
print(p9)
cat("\n\n")
print(p10)
cat("\n\n")


cat("#### Secondary pops \n\n")
  #Secondary pops
Markers_NK_6pops = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP11_6_Clusters.rds")

data.query_loop = data.query

for (i in names(Markers_NK_6pops)){
  data.query_loop = AddModuleScore(data.query_loop, features = as.list(Markers_NK_6pops[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln3 = VlnPlot(data.query_loop, features = paste0(names(Markers_NK_6pops),"1") , group.by = "predicted.id", cols= palette2, pt.size = 0)   

cat("\n\n")
print(p_Vln3)
cat("\n\n")

p11 = DotPlot(data.query_loop, features = Markers_NK_6pops$NK1A , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1A genes") + theme(axis.text = element_text(size = 20)) 
p12 = DotPlot(data.query_loop, features = Markers_NK_6pops$NK1B, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1B genes") + theme(axis.text = element_text(size = 20)) 
p13 = DotPlot(data.query_loop, features = Markers_NK_6pops$NK1C, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1C genes") + theme(axis.text = element_text(size = 20)) 
p14 = DotPlot(data.query_loop, features = Markers_NK_6pops$NKint, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NKint genes") + theme(axis.text = element_text(size = 20)) 
p15 = DotPlot(data.query_loop, features = Markers_NK_6pops$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p16 = DotPlot(data.query_loop, features = Markers_NK_6pops$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p11)
cat("\n\n")
print(p12)
cat("\n\n")
print(p13)
cat("\n\n")
print(p14)
cat("\n\n")
print(p15)
cat("\n\n")
print(p16)
cat("\n\n")



cat("### Reliability of the prediction across predicted id \n\n")

p17 = VlnPlot(data.query , features= "prediction.score.max", group.by = "predicted.id" , cols = palette2, pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p17)
cat("\n\n")



#table(data.query$predicted.id, data.query$orig.ident)
for(pop_predicted in levels(data.query$predicted.id)){
  data.query_subset = subset(data.query, subset = predicted.id == pop_predicted)
  Vlnp_loop= VlnPlot(data.query_subset , features= "prediction.score.max", group.by = "orig.ident" , pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted)
  Ridge_loop= RidgePlot(data.query_subset , features= "prediction.score.max", group.by = "orig.ident" )  & ggtitle(pop_predicted)
  
  cat("\n\n")
  print(Vlnp_loop)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop)
  cat("\n\n")
  
  Vlnp_loop2= VlnPlot(data.query_subset , features= "prediction.score.max", group.by = "Dataset" , pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted)
  Ridge_loop2= RidgePlot(data.query_subset , features= "prediction.score.max", group.by = "Dataset" )  & ggtitle(pop_predicted)
  
  cat("\n\n")
  print(Vlnp_loop2)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop2)
  cat("\n\n")
 
}


p17b = VlnPlot(data.query , features= "prediction.score.max", group.by = "orig.ident", pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black") 

cat("\n\n")
print(p17b)
cat("\n\n")

p17c = VlnPlot(data.query , features= "prediction.score.max", group.by = "predicted.id", split.by = "orig.ident", y.max = 1.1) 




cat("### Look at the spontaneous markers for the different populations predicted \n\n")

cat("#####" ,"DotPlot spontaneous markers","\n")
data.query = SetIdent(data.query, value = "predicted.id")

All_Markers = FindAllMarkers(data.query , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers
All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top10


p18 = DotPlot(data.query, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p18)
cat("\n\n")

#interactive markers table
#Make the interactive table
DEG_sample <- All_Markers
DEG_sample_filtered =(DEG_sample[DEG_sample$p_val_adj<FDRCUTOFF,])
DEG_sample_filtered <- DEG_sample_filtered[order(DEG_sample_filtered$avg_log2FC ,decreasing = T),]

cat("#####" ,"Markers Table","\n")

Nb_markers=length(unique(data.query@active.ident))
mypalette= hue_pal()(Nb_markers)
mypalette=palette2 #Vérifier si ça marche ou si ça fait bugger
print( htmltools::tagList(DT::datatable(DEG_sample_filtered,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'cluster',
    backgroundColor = DT::styleEqual(sort(unique(data.query@active.ident)),  mypalette[1:Nb_markers])
  )))


for(clust in sort(unique(Idents(data.query)))){
  cat("\n\n")
  cat("##### Cluster ", clust,  "{.tabset .tabset-fade}\n\n")
  cat("\n\n")
  DEG_clust <- DEG_sample_filtered[DEG_sample_filtered$cluster %in% clust,]
  
  for(gene in head(DEG_clust$gene)){
    cat("\n\n")
    cat("######", gene)
    cat("\n\n")
    print(VlnPlot(data.query, group.by = "predicted.id", features = gene, pt.size = 0))
    cat("\n\n")
  }
}



cat("## Look at the results of label transfer {.tabset .tabset-fade} \n\n")

#For 4 pops
cat("\n\n")
p22 = ggplot(data.query@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Query Only") + theme(text = element_text(size = 20))
print(p22)
cat("\n\n")

cat("\n\n")
p21 = ggplot(data.query@meta.data, aes(x=FirstClust, fill= predicted.id)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Query Only") + theme(text = element_text(size = 20))
print(p21)
cat("\n\n")

#For 6 pops
cat("\n\n")
p32 = ggplot(data.query2@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette) + ggtitle("Query Only") + theme(text = element_text(size = 20))
print(p32)
cat("\n\n")

cat("\n\n")
p21 = ggplot(data.query@meta.data, aes(x=FirstClust, fill= predicted.id)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Query Only") + theme(text = element_text(size = 20))
print(p21)
cat("\n\n")




# Convert to data frame for ggplot


conf_matrix <- table(data.query$FirstClust, data.query$predicted.id)

# 2. Calculate percentages
conf_matrix_percent <- sweep(conf_matrix, 1, rowSums(conf_matrix), FUN = "/") * 100

# 3. Melt the confusion matrix for ggplot
conf_matrix_melted <- reshape2::melt(conf_matrix_percent)


# 4. Plot the heatmap with percentages
plot_confu = ggplot(conf_matrix_melted, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + # Create the heatmap tiles
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) + # Reverted "RdBu"
  geom_text(aes(label=sprintf("%.1f%%", value)), color="black", size=4) + # Add text labels with percentage
  theme_minimal() + # Use a minimal theme
  labs(x = "Predicted Label", y = "True Label", fill = "Percentage") + # Label the axes and legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability

cat("\n\n")
print(plot_confu)
cat("\n\n")




conf_matrix <- table(data.query2$FirstClust, data.query2$predicted.id)

# 2. Calculate percentages
conf_matrix_percent <- sweep(conf_matrix, 1, rowSums(conf_matrix), FUN = "/") * 100

# 3. Melt the confusion matrix for ggplot
conf_matrix_melted <- reshape2::melt(conf_matrix_percent)


# 4. Plot the heatmap with percentages
plot_confu2 = ggplot(conf_matrix_melted, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + # Create the heatmap tiles
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) + # Reverted "RdBu"
  geom_text(aes(label=sprintf("%.1f%%", value)), color="black", size=4) + # Add text labels with percentage
  theme_minimal() + # Use a minimal theme
  labs(x = "Predicted Label", y = "True Label", fill = "Percentage") + # Label the axes and legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability

cat("\n\n")
print(plot_confu2)
cat("\n\n")



conf_matrix <- table(data.query2$SecondClust , data.query2$predicted.id)

# 2. Calculate percentages
conf_matrix_percent <- sweep(conf_matrix, 1, rowSums(conf_matrix), FUN = "/") * 100

# 3. Melt the confusion matrix for ggplot
conf_matrix_melted <- reshape2::melt(conf_matrix_percent)


# 4. Plot the heatmap with percentages
plot_confu2 = ggplot(conf_matrix_melted, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + # Create the heatmap tiles
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) + # Reverted "RdBu"
  geom_text(aes(label=sprintf("%.1f%%", value)), color="black", size=4) + # Add text labels with percentage
  theme_minimal() + # Use a minimal theme
  labs(x = "Predicted Label", y = "True Label", fill = "Percentage") + # Label the axes and legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability

cat("\n\n")
print(plot_confu2)
cat("\n\n")




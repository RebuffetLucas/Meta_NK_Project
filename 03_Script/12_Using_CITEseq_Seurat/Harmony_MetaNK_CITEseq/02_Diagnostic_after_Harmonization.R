## @knitr DiagnosticAfterHarmony

cat("# Diagnostic After Harmony {.tabset .tabset-fade} \n\n")

#Import object
data_all = readRDS(RDS_DIR_OUTPUT)


#Create additionnal separated objects for Tang and Meta-NK
#Focus on MetaData Only to assess the quality of batch correction
data_all_MetaNK_Only = subset(data_all, Project == "MetaNK_V2")
data_all_MetaNK_Only$orig.ident = droplevels(data_all_MetaNK_Only$orig.ident)
data_all_MetaNK_Only$Batch_Cor = droplevels(data_all_MetaNK_Only$Batch_Cor)

#Focus on Tang Only to assess the MetaNK naming across cancer
data_all_Tang_Only = subset(data_all, Project == "Tang")
data_all_Tang_Only$orig.ident = droplevels(data_all_Tang_Only$orig.ident)
data_all_Tang_Only$Batch_Cor = droplevels(data_all_Tang_Only$Batch_Cor)


#DimPlot Data viz

#Just a look at PCA before integration
p1 = DimPlot(data_all, reduction= "pca", group.by = "Project")

cat("## DimPlot \n\n")

print(p1)
cat("\n\n")

#Take a look at the result
p2 = DimPlot(data_all, group.by =  "Project")
p3 = DimPlot(data_all, label= TRUE, split.by =  "Project")
p4 = DimPlot(data_all, label= TRUE, group.by = "FirstClust", split.by = "Project")
p8 = DimPlot(data_all, label= TRUE, group.by = "SecondClust", split.by = "Project")  + NoLegend()


print(p2)
cat("\n\n")
print(p3)
cat("\n\n")
print(p4)
cat("\n\n")
print(p8)


#MetaNK alone
p5 = DimPlot(data_all_MetaNK_Only, label= TRUE)
p6 = DimPlot(data_all_MetaNK_Only, label= TRUE,  group.by = "FirstClust")
p7 = DimPlot(data_all_MetaNK_Only, label= FALSE, group.by = "orig.ident")

print(p5 + p6)
cat("\n\n")
print(  p5 + p7)
cat("\n\n")


#Barplot Data viz
  #Assess batch-correction

#Bar Diagram of Datasets, Chemistry , Orig.ident:

cat("## BarPlot \n\n")

p11 = ggplot(data_all_MetaNK_Only@meta.data, aes(x=seurat_clusters, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette) + ggtitle("MetaNK_Only")
p12 = ggplot(data_all@meta.data, aes(x=meta_histology, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + ggtitle("MetaNK and Tang")

print(  p11 +  p12)
cat("\n\n")


p13 = ggplot(data_all_MetaNK_Only@meta.data, aes(x=orig.ident, fill= seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p14 = ggplot(data_all_Tang_Only@meta.data, aes(x=Dataset, fill=seurat_clusters)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  


print(  p13 +  p14)
cat("\n\n")


#Look at markers
#Find Markers
cat("## DotPlot \n\n")

cat("\n\n")

All_Markers = FindAllMarkers(data_all , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


#saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/Integrated_Seurat_Object/DEG_12clusters.rds")
#Look at markers like Constance
markers_overCl=All_Markers
topGenes=c()
for (cl in sort(unique(Idents(data_all)))) {
  print(cl)
  markers_overCl_cl=markers_overCl[
    markers_overCl$avg_log2FC>0&
      markers_overCl$p_val_adj<=0.05&
      markers_overCl$cluster==cl,
  ]
  
  markers_overCl_cl=	
    markers_overCl_cl[order(markers_overCl_cl$p_val),]
  topGenes_cl=markers_overCl_cl$gene[
    1:min(FINDMARKERS_SHOWTOP,nrow(markers_overCl_cl))]
  topGenes[[cl]]=topGenes_cl
  print("Top genes for each clusters")
  print(	topGenes_cl)
}


topGenesb=unique(unlist(topGenes))
p10 = DotPlot(data_all, features = topGenesb , cols = "Spectral") + theme(axis.text.x = element_text(angle = 90))

cat("\n\n")

print(  p10)
cat("\n\n")

#Make the interactive table
DEG_sample <- All_Markers
DEG_sample_filtered =(DEG_sample[DEG_sample$p_val_adj<FDRCUTOFF,])
DEG_sample_filtered <- DEG_sample_filtered[order(DEG_sample_filtered$avg_log2FC ,decreasing = T),]

cat("#####" ,"Markers","\n")

Nb_markers=length(unique(data_all@active.ident))
mypalette= hue_pal()(Nb_markers)

print( htmltools::tagList(DT::datatable(DEG_sample_filtered,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'cluster',
    backgroundColor = DT::styleEqual(sort(unique(data_all@active.ident)),  mypalette[1:Nb_markers])
  )))


for(clust in sort(unique(Idents(data_all)))){
  cat("\n\n")
  cat("##### Cluster ", clust,  "{.tabset .tabset-fade}\n\n")
  cat("\n\n")
  DEG_clust <- DEG_sample_filtered[DEG_sample_filtered$cluster %in% clust,]
  
  for(gene in head(DEG_clust$gene)){
    cat("\n\n")
    cat("######", gene)
    cat("\n\n")
    print(VlnPlot(data_all, group.by = "cluster", features = gene, pt.size = 0))
    cat("\n\n")
  }
}

#Save Markers output (or not)
#saveRDS(All_Markers, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/DEG/DEG_15_clusters.rds")






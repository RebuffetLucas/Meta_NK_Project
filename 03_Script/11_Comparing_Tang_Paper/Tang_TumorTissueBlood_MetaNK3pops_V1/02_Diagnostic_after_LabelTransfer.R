## @knitr DiagnosticAfterLabelTransfer


cat("# Diagnostic After Label Transfer {.tabset .tabset-fade} \n\n")

cat("## Diag heterogeneity quality Tang Data {.tabset .tabset-fade} \n\n")


cat("### Look at sample composition before subsetting for the analysis")
#####################

#if you don't want to run all the process

#data.query = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_TumorTissueBlood_MetaNK3pops_V1/data.query_all_cells_all_tissues.rds")
#data_integrated = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_TumorTissueBlood_MetaNK3pops_V1/Ref_Integrated.rds")

data.query$meta_histology = as.factor(data.query$meta_histology)
data.query$meta_histology_relabel  =  data.query$meta_histology
levels(data.query$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(data.query$meta_histology_relabel ))




#table(PBMC_Blood_All$meta_tissue)

p02 = DimPlot(PBMC_Blood_All, group.by = "meta_tissue") 

p03 = DimPlot(PBMC_Blood_All, group.by = "meta_histology") 


cat("\n\n")
print( p02 + p03)
cat("\n\n")


table_samples = data.frame( table(PBMC_Blood_All$meta_tissue , PBMC_Blood_All$meta_histology))
colnames(table_samples) = c("Tissue" , "Cancer Type", "Number of Cells")
table_All = datatable(table_samples)

table_samples_Light = table_samples[table_samples$`Number of Cells`!=0,]



Nb_markers = length(unique(PBMC_Blood_All$meta_histology))
mypalette= hue_pal()(Nb_markers)

cat("\n\n")

print( htmltools::tagList(DT::datatable(table_samples_Light,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'Tissue',
    backgroundColor = DT::styleEqual(sort(unique(PBMC_Blood_All$meta_histology)),  mypalette[1:Nb_markers])
  )))


cat("\n\n")
#####################


cat("### Tang and Meta NK data raw \n\n")

#####################

  #Global UMAP
p01 = DimPlot(PBMC_Meta_NK_V2, cols = palette3) 

p02 = DimPlot(PBMC_Blood_All, group.by = "meta_tissue") 

p03 = DimPlot(PBMC_Blood_All, group.by = "meta_histology") 



cat("\n\n")
print(p01)
cat("\n\n")

cat("\n\n")
print( p02 + p03)
cat("\n\n")
#####################


cat("### Diag Row Tang data Counts and Features \n\n")
#Diag of Counts an Features
#####################
PBMC_Blood_All$nCount = colSums(x = PBMC_Blood_All, slot = "counts")  # nCount_RNA
PBMC_Blood_All$nFeature = colSums(x = GetAssayData(object = PBMC_Blood_All, slot = "counts") > 0)  # nFeatureRNA
PBMC_Blood_All$tissue_Firstclust = paste0(PBMC_Blood_All$meta_tissue , "_", PBMC_Blood_All$celltype )

p04 = VlnPlot(PBMC_Blood_All, features= "nCount", group.by = "meta_histology")
p05 = VlnPlot(PBMC_Blood_All, features= "nCount", group.by = "datasets")

p06 = VlnPlot(PBMC_Blood_All, features= "nCount", group.by = "tissue_Firstclust")

cat("\n\n")
print(p04)
cat("\n\n")

cat("\n\n")
print(p05)
cat("\n\n")

cat("\n\n")
print(p06)
cat("\n\n")

#####################


cat("### Look at sample after subsetting for the analysis")
#####################

#UMAP data viz
Cells_DataQuery = Cells(data.query)
Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)


p02 = DimPlot(PBMC_Blood_All, group.by = "meta_tissue", cells =Cells_DataQuery ) 
p03 = DimPlot(PBMC_Blood_All, group.by = "meta_histology", cells = Cells_DataQuery) 


cat("\n\n")
print( p02 + p03)
cat("\n\n")



#table(PBMC_Blood_All$meta_tissue)

table_samples = data.frame( table(data.query$meta_tissue , data.query$meta_histology))
colnames(table_samples) = c("Tissue" , "Cancer Type", "Number of Cells")
table_All2 = datatable(table_samples)

table_samples_Light = table_samples[table_samples$`Number of Cells`!=0,]




Nb_markers = length(unique(data.query$meta_histology))
mypalette= hue_pal()(Nb_markers)


cat("\n\n")

print( htmltools::tagList(DT::datatable(table_samples_Light,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'Tissue',
    backgroundColor = DT::styleEqual(sort(unique(data.query$meta_histology)),  mypalette[1:Nb_markers])
  )))


cat("\n\n")

#####################

cat("## Vizualization of the reference {.tabset .tabset-fade} \n\n")
#####################

  #Have a look at the reference built with the 4 Dataset of Meta-NK

p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "Dataset")
p2 <- DimPlot(data_integrated, reduction = "umap", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette3) + NoLegend()

cat("\n\n")
print( p1 + p2 )

cat("\n\n")

p4 <- DimPlot(data_integrated, reduction = "umap", split.by = "Dataset", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette3) + NoLegend()

print(p4)
cat("\n\n")


cat("## Assessing quality of Label Transfer {.tabset .tabset-fade} \n\n")

cat("### Look at the common markers for the different populations {.tabset .tabset-fade} \n\n")

data.query$predicted.id = as.factor(data.query$predicted.id)

data.query = SetIdent(data.query, value = "predicted.id")

#####################

cat("#### Main pops based on CITEseq \n\n")

#####################

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
data.query2 = data.query
data.query2 = AddModuleScore(data.query2, features = list(Top_NK1$gene), pool= NULL ,name= "NK1" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK2$gene), pool= NULL ,name= "NK2" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK3$gene), pool= NULL ,name= "NK3" , seed=19)

list_data.query_tissue = SplitObject(data.query2 , split.by = "meta_tissue")

for (data_query_tissue in list_data.query_tissue){
 p_Vln2 = VlnPlot(data_query_tissue, features = c("NK11", "NK21", "NK31") , group.by = "predicted.id", cols= palette3, pt.size = 0)   & stat_summary(fun.data=data_summary,color="black" ) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))
 
   
cat("\n\n")
print(p_Vln2)
cat("\n\n")
  
p8 = DotPlot(data_query_tissue, features = Top_NK1$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) + ggtitle( unique(data_query_tissue$meta_tissue) )
p9 = DotPlot(data_query_tissue, features = Top_NK2$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) + ggtitle( unique(data_query_tissue$meta_tissue) )
p10 = DotPlot(data_query_tissue, features = Top_NK3$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) + ggtitle( unique(data_query_tissue$meta_tissue) )

cat("\n\n")
print(p8)
cat("\n\n")
print(p9)
cat("\n\n")
print(p10)
cat("\n\n")



}
#####################





cat("#### Secondary pops \n\n")
#####################

  #Secondary pops
Markers_NK_6pops = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP11_6_Clusters.rds")

data.query2 = data.query

for (i in names(Markers_NK_6pops)){
  data.query2 = AddModuleScore(data.query2, features = as.list(Markers_NK_6pops[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln3 = VlnPlot(data.query2, features = paste0(names(Markers_NK_6pops),"1") , group.by = "predicted.id", cols= palette3, pt.size = 0)   

cat("\n\n")
print(p_Vln3)
cat("\n\n")

p11 = DotPlot(data.query, features = Markers_NK_6pops$NK1A , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1A genes") + theme(axis.text = element_text(size = 20)) 
p12 = DotPlot(data.query, features = Markers_NK_6pops$NK1B, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1B genes") + theme(axis.text = element_text(size = 20)) 
p13 = DotPlot(data.query, features = Markers_NK_6pops$NK1C, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1C genes") + theme(axis.text = element_text(size = 20)) 
p14 = DotPlot(data.query, features = Markers_NK_6pops$NKint, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NKint genes") + theme(axis.text = element_text(size = 20)) 
p15 = DotPlot(data.query, features = Markers_NK_6pops$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p16 = DotPlot(data.query, features = Markers_NK_6pops$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

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
#####################


cat("### Reliability of the prediction across predicted id \n\n")


#####################

p17 = VlnPlot(data.query , features= "prediction.score.max", group.by = "predicted.id" , cols = palette3, pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black") & ggtitle(" All Data Included")

list_data.query_tissue = SplitObject(data.query , split.by = "meta_tissue")

cat("\n\n")
print(p17)
cat("\n\n")

for (data_query_tissue in list_data.query_tissue){
  p17_tissue = VlnPlot(data_query_tissue , features= "prediction.score.max", group.by = "predicted.id" , cols = palette3, pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black")& plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))
  cat("\n\n")
  print(p17_tissue)
  cat("\n\n")
  
}



#table(data.query$predicted.id, data.query$meta_tissue)
for (data_query_tissue in list_data.query_tissue){
for(pop_predicted in levels(data_query_tissue$predicted.id)){
  data.query_subset = subset(data_query_tissue, subset = predicted.id == pop_predicted)
  
  data.query_subset@meta.data = droplevels(data.query_subset@meta.data)
  
  Vlnp_loop3= VlnPlot(data.query_subset , features= "prediction.score.max", group.by = "meta_histology" , pt.size = 0) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))
  Ridge_loop3= RidgePlot(data.query_subset , features= "prediction.score.max", sort =TRUE ,group.by = "meta_histology" )  & ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))
  
  cat("\n\n")
  print(Vlnp_loop3)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop3)
  cat("\n\n")
}

}


for (data_query_tissue in list_data.query_tissue){

p17b = VlnPlot(data_query_tissue , features= "prediction.score.max", group.by = "meta_histology", pt.size = 0) & stat_summary(fun.data=data_summary,color="black")  & ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))


cat("\n\n")
print(p17b)
cat("\n\n")

p17c = VlnPlot(data_query_tissue , features= "prediction.score.max", group.by = "predicted.id", split.by = "meta_histology") & ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))


cat("\n\n")
print(p17c)
cat("\n\n")

p17d = VlnPlot(data_query_tissue , features= "prediction.score.max", group.by = "meta_histology", split.by = "predicted.id" , cols = palette3) & ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15))

p17d

cat("\n\n")
print(p17d)
cat("\n\n")

}



p17b = VlnPlot(data.query , features= "prediction.score.max", group.by = "meta_histology", pt.size = 0) & stat_summary(fun.data=data_summary,color="black") 

cat("\n\n")
print(p17b)
cat("\n\n")

p17c = VlnPlot(data.query , features= "prediction.score.max", group.by = "predicted.id", split.by = "meta_histology") 

cat("\n\n")
print(p17c)
cat("\n\n")

p17d = VlnPlot(data.query , features= "prediction.score.max", group.by = "meta_histology", split.by = "predicted.id" , cols = palette3)


cat("\n\n")
print(p17d)
cat("\n\n")
#####################



cat("### Look at the spontaneous markers for the different populations predicted \n\n")

cat("#####" ,"DotPlot spontaneous markers","\n")

#####################


data.query = SetIdent(data.query, value = "predicted.id")


for (data_query_tissue in list_data.query_tissue){
  
  
  
  All_Markers = FindAllMarkers(data_query_tissue , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)
  
  #Look at markers
  All_Markers %>%
    group_by(cluster) %>%
    filter(!grepl("RPL|RPS|MT-", gene)) %>%
    top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top10
  
  
  p18_loop = DotPlot(data.query, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(size = 20)) + ggtitle( unique(data_query_tissue$meta_tissue) )
  
  cat("\n\n")
  print(p18_loop)
  cat("\n\n")
  
  
}

All_Markers = FindAllMarkers(data.query , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


#Look at markers
All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top10


p18 = DotPlot(data.query, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(size = 20)) + ggtitle("With all data included")

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
mypalette=palette3 #Vérifier si ça marche ou si ça fait bugger
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

#####################


cat("## Look at the results of label transfer {.tabset .tabset-fade} \n\n")

#####################

#Re-inject UMAP coordinates in data.query

Cells_DataQuery = Cells(data.query)
Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
PBMC_Blood_All_subs_data_query$predicted.id = data.query$predicted.id
PBMC_Blood_All_subs_data_query$meta_histology_relabel  = PBMC_Blood_All_subs_data_query$meta_histology
levels(PBMC_Blood_All_subs_data_query$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(PBMC_Blood_All_subs_data_query$meta_histology_relabel ))

#PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All_subs_data_query, subset = meta_tissue == "Tumor") #To look only at the tumor
title = paste(unique(PBMC_Blood_All_subs_data_query$meta_tissue), collapse = " ")



cat("### All data Tang ", title)


Clusters_Proportions = prop.table( table( PBMC_Blood_All_subs_data_query$meta_histology , PBMC_Blood_All_subs_data_query$predicted.id) , margin = 1)
df  = as.matrix(Clusters_Proportions)
df = as.data.frame.matrix(Clusters_Proportions)
df$meta_histology <- rownames(df)
df_filtered <- df[!is.na(df$NK2), ]
df_ordered <- df_filtered[order(df_filtered$NK2), ]


ORDER_CANCER = rownames(df_ordered)
ORDER_CANCER_RELABEL = sub(".*\\((.*)\\).*", "\\1", ORDER_CANCER )

PBMC_Blood_All_subs_data_query$meta_histology = droplevels(PBMC_Blood_All_subs_data_query$meta_histology)
PBMC_Blood_All_subs_data_query$meta_histology  = factor(PBMC_Blood_All_subs_data_query$meta_histology  , levels= ORDER_CANCER)

PBMC_Blood_All_subs_data_query$meta_histology_relabel = droplevels(PBMC_Blood_All_subs_data_query$meta_histology_relabel )
PBMC_Blood_All_subs_data_query$meta_histology_relabel  = factor(PBMC_Blood_All_subs_data_query$meta_histology_relabel , levels= ORDER_CANCER_RELABEL )




  cat("\n\n")

p22 = ggplot(PBMC_Blood_All_subs_data_query@meta.data, aes(x=celltype, fill= predicted.id)) + geom_bar(position="fill")  + scale_fill_manual(values= palette3)+ ggtitle(unlist(paste0("Tang Only", title)) ) +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 
  
                                            
cat("\n\n")
print(p22)
cat("\n\n")


p21 = ggplot(PBMC_Blood_All_subs_data_query@meta.data, aes(x=meta_histology, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
  ggtitle(unlist(paste0("Tang Only", title)) )  +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 


p21bis = ggplot(PBMC_Blood_All_subs_data_query@meta.data, aes(x=meta_histology_relabel, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
  ggtitle(unlist(paste0("Tang Only", title)) )  +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) +
  scale_x_discrete(labels = label_wrap(20)) 

cat("\n\n")
print(p21)
cat("\n\n")


cat("\n\n")
print(p21bis)
cat("\n\n")




cat("\n\n")
p_pred_id = DimPlot(PBMC_Blood_All_subs_data_query , group.by = "predicted.id", cols= palette3)  + ggtitle("Tang Only") + theme(text = element_text(size = 10)) + NoAxes()
print(p_pred_id)
cat("\n\n")






#####################


#Make a list
list_data.query_tissue = SplitObject(data.query , split.by = "meta_tissue")

#####################

for (data_query_tissue in list_data.query_tissue){
  cat("\n\n")
  cat("### Tang Only" , unique(data_query_tissue$meta_tissue)   )
  title = paste(unique(data_query_tissue$meta_tissue), collapse = " ")
  
  
  
  cat("\n\n")

  
  p22 = ggplot(data_query_tissue@meta.data, aes(x= SecondClust, fill= predicted.id)) + geom_bar(position="fill")  + scale_fill_manual(values= palette3)+ ggtitle(unlist(paste0("Tang Only", title)) ) +
    theme_classic()+
    theme(text = element_text(size = 20), axis.text.x = element_text(hjust=1, angle = 45) ) +
    scale_x_discrete(labels = label_wrap(20)) 
  
  cat("\n\n")
  print(p22)
  cat("\n\n")
  
  
  p21 = ggplot(data_query_tissue@meta.data, aes(x=meta_histology, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
    ggtitle(unlist(paste0("Tang Only", title)) )  +
    theme_classic()+
    theme(text = element_text(size = 20), axis.text.x = element_text(hjust=1, angle = 45) ) +
    scale_x_discrete(labels = label_wrap(20)) 
  
  
  p21bis = ggplot(data_query_tissue@meta.data, aes(x=meta_histology_relabel, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
    ggtitle(unlist(paste0("Tang Only", title)) )  +
    theme_classic()+
    theme(text = element_text(size = 20), axis.text.x = element_text(hjust=1, angle = 45) ) 
    
  
  cat("\n\n")
  print(p21)
  cat("\n\n")
  

  
  cat("\n\n")
  print(p21bis)
  cat("\n\n")
  
  
cat("\n\n")
p_pred_id = DimPlot(data_query_tissue , group.by = "predicted.id", cols= palette3)  + ggtitle("Tang Only") + theme(text = element_text(size = 20))
print(p_pred_id)
cat("\n\n")


 
  
  
  
  #Make a UMAP of the predicted id
  Cells_DataQuery = Cells(data_query_tissue)
  Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
  PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
  PBMC_Blood_All_subs_data_query$predicted.id = data_query_tissue$predicted.id
  
  cat("\n\n")
  p_pred_id = DimPlot(PBMC_Blood_All_subs_data_query , group.by = "predicted.id", cols= palette3)  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", unique(data_query_tissue$meta_tissue), " predicted id")) + theme(text = element_text(size = 20))
  print(p_pred_id)
  cat("\n\n")
  
}

cat("### Meta NK as a comparison")

cat("\n\n")
p23 = ggplot(PBMC_Meta_NK_V2@meta.data, aes(x=orig.ident, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + scale_fill_manual(values= palette3) + ggtitle("Meta NK per Sample") + theme(text = element_text(size = 20))

print(p23)
cat("\n\n")

cat("\n\n")
p24 = ggplot(PBMC_Meta_NK_V2@meta.data, aes(x=Dataset, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + scale_fill_manual(values= palette3)  + ggtitle("Meta NK per Dataset") + theme(text = element_text(size = 20))

print(p24)
cat("\n\n")
#####################


cat("## Scoring on Tang UMAP {.tabset .tabset-fade} \n\n")
#####################

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
data.query2 = data.query
data.query2 = AddModuleScore(data.query2, features = list(Top_NK1$gene), pool= NULL ,name= "NK1" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK2$gene), pool= NULL ,name= "NK2" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK3$gene), pool= NULL ,name= "NK3" , seed=19)

list_data.query_tissue = SplitObject(data.query2 , split.by = "meta_tissue")



cat("### All data Tang")

Cells_DataQuery = Cells(data.query2)
Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
PBMC_Blood_All_subs_data_query$NK1 = data.query2$NK11
PBMC_Blood_All_subs_data_query$NK2 = data.query2$NK21
PBMC_Blood_All_subs_data_query$NK3 = data.query2$NK31

cat("\n\n")
p_score = FeaturePlot(PBMC_Blood_All_subs_data_query , feature = "NK1")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", "All data ", " NK1 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p_score)
cat("\n\n")

cat("\n\n")
p_score = FeaturePlot(PBMC_Blood_All_subs_data_query , feature = "NK2")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", "All data ", " NK2 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p_score)
cat("\n\n")

cat("\n\n")
p_score = FeaturePlot(PBMC_Blood_All_subs_data_query , feature = "NK3")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", "All data ", " NK3 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p_score)
cat("\n\n")


for (data_query_tissue in list_data.query_tissue){
  
  Cells_DataQuery = Cells(data_query_tissue)
  Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
  PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
  PBMC_Blood_All_subs_data_query$NK1 = data_query_tissue$NK11
  PBMC_Blood_All_subs_data_query$NK2 = data_query_tissue$NK21
  PBMC_Blood_All_subs_data_query$NK3 = data_query_tissue$NK31
  
  cat("\n\n")
  cat("### All data" , unique(data_query_tissue$meta_tissue)   )
  
  cat("\n\n")
  p_score = FeaturePlot(PBMC_Blood_All_subs_data_query , feature = "NK1")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", unique(PBMC_Blood_All_subs_data_query$meta_tissue), " NK1 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  cat("\n\n")
  p_score = FeaturePlot(PBMC_Blood_All_subs_data_query , feature = "NK2")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", unique(PBMC_Blood_All_subs_data_query$meta_tissue), " NK2 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  cat("\n\n")
  p_score = FeaturePlot(PBMC_Blood_All_subs_data_query , feature = "NK3")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", unique(PBMC_Blood_All_subs_data_query$meta_tissue), " NK3 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
}
#####################


cat("## Look at the results on clusters on tissue and cancer type one by one {.tabset .tabset-fade} \n\n")

#####################

Cells_DataQuery = Cells(data.query2)
Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
PBMC_Blood_All_subs_data_query@meta.data = data.query2@meta.data

rownames(PBMC_Blood_All_subs_data_query@meta.data) = str_sub(rownames(PBMC_Blood_All_subs_data_query@meta.data) , 6,-1)


list_data.query_Tissue_Histo = SplitObject(PBMC_Blood_All_subs_data_query , split.by = "tissue_Histo")


for (i in 1:length(list_data.query_Tissue_Histo)) {
  Name_Sample = names(list_data.query_Tissue_Histo[i])
  Number_Cells_Sample = length(Cells(list_data.query_Tissue_Histo[[i]]))
  
  cat("###  " , Name_Sample   )
  
  cat("\n\n")
  print(paste0("Number of cells:   " , Number_Cells_Sample))
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "predicted.id") + ggtitle(paste0(Name_Sample , " predicted in ", " Tang UMAP")))
  cat("\n\n")
  
  list_data.query_Tissue_Histo[[i]] <- FindVariableFeatures(list_data.query_Tissue_Histo[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  
  list_data.query_Tissue_Histo[[i]] <- ScaleData(list_data.query_Tissue_Histo[[i]] , verbose = FALSE)
  list_data.query_Tissue_Histo[[i]] <- RunPCA(list_data.query_Tissue_Histo[[i]], npcs = 30, verbose = FALSE)
  list_data.query_Tissue_Histo[[i]] <- RunUMAP(list_data.query_Tissue_Histo[[i]], reduction = "pca", dims = 1:30)
  
  #Global UMAP with predicted annotation
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "predicted.id", cols= palette3) + ggtitle(paste0(Name_Sample , " predicted in ", " New UMAP")))
  cat("\n\n")
  
  
  #Global UMAP scoring with NK1,2 and 3 signatures
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK11")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", Name_Sample, " NK1 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK21")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", Name_Sample, " NK2 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK31")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", Name_Sample, " NK3 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
}

#####################

cat("## Look at the results of prediction and scores on tissue one by one ( no correction) {.tabset .tabset-fade} \n\n")

#####################


Cells_DataQuery = Cells(data.query2)
Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
PBMC_Blood_All_subs_data_query@meta.data = data.query2@meta.data

rownames(PBMC_Blood_All_subs_data_query@meta.data) = str_sub(rownames(PBMC_Blood_All_subs_data_query@meta.data) , 6,-1)


list_data.query_Tissue_Histo = SplitObject(PBMC_Blood_All_subs_data_query , split.by = "meta_tissue")


for (i in 1:length(list_data.query_Tissue_Histo)) {
  Name_Sample = names(list_data.query_Tissue_Histo[i])
  Number_Cells_Sample = length(Cells(list_data.query_Tissue_Histo[[i]]))
  
  cat("###  " , Name_Sample   )
  
  cat("\n\n")
  print(paste0("Number of cells:   " , Number_Cells_Sample))
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "predicted.id", cols= palette3) + ggtitle(paste0(Name_Sample , " predicted in ", " Tang UMAP")))
  cat("\n\n")
  
  list_data.query_Tissue_Histo[[i]] <- FindVariableFeatures(list_data.query_Tissue_Histo[[i]], selection.method = "vst", 
                                                            nfeatures = 2000, verbose = FALSE)
  
  list_data.query_Tissue_Histo[[i]] <- ScaleData(list_data.query_Tissue_Histo[[i]] , verbose = FALSE)
  list_data.query_Tissue_Histo[[i]] <- RunPCA(list_data.query_Tissue_Histo[[i]], npcs = 30, verbose = FALSE)
  list_data.query_Tissue_Histo[[i]] <- RunUMAP(list_data.query_Tissue_Histo[[i]], reduction = "pca", dims = 1:30)
  
  #Global UMAP with predicted annotation
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "predicted.id", cols= palette3) + ggtitle(paste0(Name_Sample , " predicted in ", " New UMAP")))
  cat("\n\n")
  
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "meta_histology") + ggtitle(paste0(Name_Sample , " predicted in ", " New UMAP")))
  cat("\n\n")
  
  #Global UMAP scoring with NK1,2 and 3 signatures
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK11")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", Name_Sample, " NK1 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK21")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", Name_Sample, " NK2 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK31")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only", Name_Sample, " NK3 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
}
#####################


cat("## Look at the results of prediction and scores on tissue one by one (with harmony correction) {.tabset .tabset-fade} \n\n")

#####################

Cells_DataQuery = Cells(data.query2)
Cells_DataQuery = str_sub(Cells_DataQuery , 6,-1)
PBMC_Blood_All_subs_data_query = subset(PBMC_Blood_All, cells= Cells_DataQuery)
PBMC_Blood_All_subs_data_query@meta.data = data.query2@meta.data

rownames(PBMC_Blood_All_subs_data_query@meta.data) = str_sub(rownames(PBMC_Blood_All_subs_data_query@meta.data) , 6,-1)


list_data.query_Tissue_Histo = SplitObject(PBMC_Blood_All_subs_data_query , split.by = "meta_tissue")


for (i in 1:length(list_data.query_Tissue_Histo)) {
  Name_Sample = names(list_data.query_Tissue_Histo[i])
  Number_Cells_Sample = length(Cells(list_data.query_Tissue_Histo[[i]]))
  
  cat("###  " , Name_Sample   )
  
  cat("\n\n")
  print(paste0("Number of cells:   " , Number_Cells_Sample))
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "predicted.id", cols= palette3) + ggtitle(paste0(Name_Sample , " predicted in ", " Tang UMAP")))
  cat("\n\n")
  
  list_data.query_Tissue_Histo[[i]] <- FindVariableFeatures(list_data.query_Tissue_Histo[[i]], selection.method = "vst", 
                                                            nfeatures = 2000, verbose = FALSE)
  
  list_data.query_Tissue_Histo[[i]] <- ScaleData(list_data.query_Tissue_Histo[[i]] , verbose = FALSE)
  list_data.query_Tissue_Histo[[i]] <- RunPCA(list_data.query_Tissue_Histo[[i]], npcs = 30, verbose = FALSE)
  list_data.query_Tissue_Histo[[i]] = list_data.query_Tissue_Histo[[i]] %>%   RunHarmony(c("meta_histology"), plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:30 ) 
  list_data.query_Tissue_Histo[[i]] <- RunUMAP(list_data.query_Tissue_Histo[[i]], reduction = "harmony", dims = 1:30)
  
  #Global UMAP with predicted annotation
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "predicted.id", cols= palette3) + ggtitle(paste0(Name_Sample , " predicted in ", " New UMAP")))
  cat("\n\n")
  
  
  cat("\n\n")
  print(DimPlot(list_data.query_Tissue_Histo[[i]], group.by = "meta_histology") + ggtitle(paste0(Name_Sample , " predicted in ", " New UMAP")))
  cat("\n\n")
  
  #Global UMAP scoring with NK1,2 and 3 signatures
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK11")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only ", Name_Sample, " NK1 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK21")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only ", Name_Sample, " NK2 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
  
  cat("\n\n")
  p_score = FeaturePlot(list_data.query_Tissue_Histo[[i]] , feature = "NK31")  + theme(text = element_text(size = 20)) + ggtitle(paste0("Tang Only ", Name_Sample, " NK3 score")) + theme(text = element_text(size = 20)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p_score)
  cat("\n\n")
}


#####################











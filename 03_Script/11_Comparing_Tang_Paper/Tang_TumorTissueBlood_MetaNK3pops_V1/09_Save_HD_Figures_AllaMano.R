#HD Plots à  la mano parce que ça m'a saoûlé



#Plot for quality of Data prediction

#Barplot proportions

cat("## Look at the results of label transfer {.tabset .tabset-fade} \n\n")

#####################

#Load the Data
data.query = readRDS( DATA_QUERY_SAVE_PATH)
data_integrated = readRDS(DATA_REF_SAVE_PATH)



#Prepare the objects
data.query$meta_histology= as.factor(data.query$meta_histology)
data.query$meta_histology_relabel  = data.query$meta_histology
levels(data.query$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(data.query$meta_histology_relabel ))

data.query$predicted.id =  as.factor(data.query$predicted.id)
data.query$SecondClust  = as.factor(data.query$SecondClust)


list_data.query_tissue = SplitObject(data.query , split.by = "meta_tissue")


for (data_query_tissue in list_data.query_tissue){

title = paste(unique(data_query_tissue$meta_tissue), collapse = " ")



cat("### All data Tang ", title)

Clusters_Proportions = prop.table( table( data_query_tissue$meta_histology , data_query_tissue$predicted.id) , margin = 1)
df  = as.matrix(Clusters_Proportions)
df = as.data.frame.matrix(Clusters_Proportions)
df$meta_histology <- rownames(df)
df_filtered <- df[!is.na(df$NK2), ]
df_ordered <- df_filtered[order(df_filtered$NK2), ]


ORDER_CANCER = rownames(df_ordered)
ORDER_CANCER_RELABEL = sub(".*\\((.*)\\).*", "\\1", ORDER_CANCER )

data_query_tissue$meta_histology = droplevels(data_query_tissue$meta_histology)
data_query_tissue$meta_histology  = factor(data_query_tissue$meta_histology  , levels= ORDER_CANCER)

data_query_tissue$meta_histology_relabel = droplevels(data_query_tissue$meta_histology_relabel )
data_query_tissue$meta_histology_relabel  = factor(data_query_tissue$meta_histology_relabel , levels= ORDER_CANCER_RELABEL )

if (title=="Tumor"){
saveRDS(ORDER_CANCER, file=paste0(PATH_ANALYSIS_OUTPUT, "/order_cancer_Tumor.rds"))
saveRDS(ORDER_CANCER_RELABEL, file=paste0(PATH_ANALYSIS_OUTPUT, "/order_cancer_Tumor_relabeled.rds"))
}

cat("\n\n")

p22 = ggplot(data_query_tissue@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill")  + scale_fill_manual(values= palette3)+ ggtitle(unlist(paste0("Tang Only", title)) ) +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 


cat("\n\n")
print(p22)
cat("\n\n")


p21 = ggplot(data_query_tissue@meta.data, aes(x=meta_histology, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
  ggtitle(unlist(paste0("Tang Only", title)) )  +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 



p21bis = ggplot(data_query_tissue@meta.data, aes(x=meta_histology_relabel, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
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


png(file=paste0(FILE_OUTPUT_PLOTS , "/MetaNKvsTangClusters",title,".png"), width = 15, height = 10,  units = "cm", res=600 )
print(p22)
dev.off()


png(file=paste0(FILE_OUTPUT_PLOTS , "/Proportions_PerCancer",title,".png"), width = 15, height = 10,  units = "cm", res=600 )
print(p21)
dev.off()


png(file=paste0(FILE_OUTPUT_PLOTS , "/Proportions_PerCancer_Abreviated",title,".png"), width = 15, height = 10,  units = "cm", res=600 )
print(p21bis)
dev.off()



Clusters_Proportions = prop.table( table( data_query_tissue$SecondClust , data_query_tissue$predicted.id) , margin = 1)
df  = as.matrix(Clusters_Proportions)
df = as.data.frame.matrix(Clusters_Proportions)
df$SecondClust <- rownames(df)
df_filtered <- df[!is.na(df$NK2), ]
df_ordered <- df_filtered[order(df_filtered$NK2), ]


ORDER_CANCER = rownames(df_ordered)
ORDER_CANCER_RELABEL = sub(".*\\((.*)\\).*", "\\1", ORDER_CANCER )


data_query_tissue$SecondClust = droplevels(data_query_tissue$SecondClust)
data_query_tissue$SecondClust  = factor(data_query_tissue$SecondClust  , levels= ORDER_CANCER)



cat("\n\n")

p25 = ggplot(data_query_tissue@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill")  + scale_fill_manual(values= palette3)+ ggtitle(unlist(paste0("Tang Only", title)) ) +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 

png(file=paste0(FILE_OUTPUT_PLOTS , "/MetaNKvsTangClustersOrdered",title,".png"), width = 15, height = 10,  units = "cm", res=600 )
print(p25)
dev.off()


}



#Plot for quality of Data prediction


cat("### Reliability of the prediction across predicted id \n\n")


#####################

p17 = VlnPlot(data.query , features= "prediction.score.max", group.by = "predicted.id" , cols = palette3, pt.size = 0, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black") & ggtitle(" All Data Included")

data.query$meta_tissue = as.factor(data.query$meta_tissue)
data.query$predicted.id = as.factor(data.query$predicted.id)
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
    
    Vlnp_loop3= VlnPlot(data.query_subset , features= "prediction.score.max", group.by = "meta_histology_relabel" , pt.size = 0, sort = "decreasing") & 
      stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) &
      theme(plot.title = element_text(hjust=0.5, size= 6))   &      theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45, size = 7), axis.text.y = element_text(hjust=1,  size = 7) ) 
    
    FILE_SAVE= paste0(FILE_OUTPUT_PLOTS,"/VlnPlot_prediction_score",as.character(unique(data_query_tissue$meta_tissue)), as.character(pop_predicted) ,".png")
    
    png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
    print(Vlnp_loop3)
    dev.off()
    
    
    
    Ridge_loop3= RidgePlot(data.query_subset , features= "prediction.score.max", sort = "decreasing" ,group.by = "meta_histology_relabel" ) &
     ggtitle(pop_predicted) & plot_annotation(title = unique(data_query_tissue$meta_tissue)) & theme(plot.title = element_text(hjust=0.5, size= 15)) & 
      theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45, size = 7), axis.text.y = element_text(hjust=1,  size = 7))
    
    FILE_SAVE= paste0(FILE_OUTPUT_PLOTS,"/RidgePlot_prediction_score",as.character(unique(data_query_tissue$meta_tissue)), as.character(pop_predicted) ,".png")
    
    png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
    print(Ridge_loop3)
    dev.off()
    
    
    cat("\n\n")
    print(Vlnp_loop3)
    cat("\n\n")
    
    cat("\n\n")
    print(Ridge_loop3)
    cat("\n\n")
  }
  
}





p17b = VlnPlot(data.query , features= "prediction.score.max", group.by = "meta_histology_relabel", pt.size = 0, sort = "decreasing" ) & stat_summary(fun.data=data_summary,color="black") & ggtitle( "All NK") & plot_annotation(title = paste(unique(data.query$meta_tissue), collapse = " ") ) &
  theme(plot.title = element_text(hjust=0.5, size= 6))   &      theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45, size = 7), axis.text.y = element_text(hjust=1,  size = 7) )

cat("\n\n")
print(p17b)
cat("\n\n")



png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p17b)
dev.off()


p17c = VlnPlot(list_data.query_tissue$Tumor , features= "prediction.score.max", group.by = "meta_histology_relabel", pt.size = 0, sort = "decreasing" ) & stat_summary(fun.data=data_summary,color="black") & ggtitle( "All NK") & plot_annotation(title = paste(unique(list_data.query_tissue$Tumor$meta_tissue), collapse = " ") ) &
  theme(plot.title = element_text(hjust=0.5, size= 6))   &      theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45, size = 7), axis.text.y = element_text(hjust=1,  size = 7) )

FILE_SAVE= paste0(FILE_OUTPUT_PLOTS,"/VlnPlot_prediction_score_AllNK_Tumor" ,".png")

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p17c)
dev.off()


p17d = VlnPlot(list_data.query_tissue$Normal , features= "prediction.score.max", group.by = "meta_histology_relabel", pt.size = 0, sort = "decreasing" ) & stat_summary(fun.data=data_summary,color="black") & ggtitle( "All NK") & plot_annotation(title = paste(unique(list_data.query_tissue$Normal$meta_tissue), collapse = " ") ) &
  theme(plot.title = element_text(hjust=0.5, size= 6))   &      theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45, size = 7), axis.text.y = element_text(hjust=1,  size = 7) )

FILE_SAVE= paste0(FILE_OUTPUT_PLOTS,"/VlnPlot_prediction_score_AllNK_Normal" ,".png")

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p17d)
dev.off()




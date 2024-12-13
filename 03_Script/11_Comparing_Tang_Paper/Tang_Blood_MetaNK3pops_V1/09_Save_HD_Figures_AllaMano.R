#HD Plots à  la mano parce que ça m'a saoûlé


#Barplot proportions

cat("## Look at the results of label transfer {.tabset .tabset-fade} \n\n")

#####################

#Load the Data
data.query = readRDS(DATA_QUERY_SAVE_PATH)
data_integrated = readRDS(DATA_REF_SAVE_PATH)



#Re-inject UMAP coordinates in data.query
data.query$meta_histology= as.factor(data.query$meta_histology)
data.query$meta_histology_relabel  = data.query$meta_histology
levels(data.query$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(data.query$meta_histology_relabel ))

data.query$predicted.id =  as.factor(data.query$predicted.id)




title = paste(unique(data.query$meta_tissue), collapse = " ")



cat("### All data Tang ", title)


Clusters_Proportions = prop.table( table( data.query$meta_histology , data.query$predicted.id) , margin = 1)
df  = as.matrix(Clusters_Proportions)
df = as.data.frame.matrix(Clusters_Proportions)
df$meta_histology <- rownames(df)
df_filtered <- df[!is.na(df$NK2), ]
df_ordered <- df_filtered[order(df_filtered$NK2), ]


ORDER_CANCER = rownames(df_ordered)
ORDER_CANCER_RELABEL = sub(".*\\((.*)\\).*", "\\1", ORDER_CANCER )

data.query$meta_histology = droplevels(data.query$meta_histology)
data.query$meta_histology  = factor(data.query$meta_histology  , levels= ORDER_CANCER)

data.query$meta_histology_relabel = droplevels(data.query$meta_histology_relabel )
data.query$meta_histology_relabel  = factor(data.query$meta_histology_relabel , levels= ORDER_CANCER_RELABEL )




cat("\n\n")

p22 = ggplot(data.query@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill")  + scale_fill_manual(values= palette3)+ ggtitle(unlist(paste0("Tang Only", title)) ) +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 


cat("\n\n")
print(p22)
cat("\n\n")


p21 = ggplot(data.query@meta.data, aes(x=meta_histology, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
  ggtitle(unlist(paste0("Tang Only", title)) )  +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 



p21bis = ggplot(data.query@meta.data, aes(x=meta_histology_relabel, fill= predicted.id)) + geom_bar(position="fill")   + scale_fill_manual(values= palette3) + 
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


#Saving HD Figures
png(file=paste0(FILE_OUTPUT_PLOTS,"/Barplot_TangBlood_proportions_TangnamingvsMetaNK.png"), width = 15, height = 10,  units = "cm", res=600 )
p22
dev.off()


png(file=paste0(FILE_OUTPUT_PLOTS,"/Barplot_TangBlood_proportionspercancer.png"), width = 15, height = 10,  units = "cm", res=600 )
p21
dev.off()


png(file=paste0(FILE_OUTPUT_PLOTS,"/Barplot_TangBlood_proportionspercancer_abreviation.png"), width = 15, height = 10,  units = "cm", res=600 )
p21bis
dev.off()


#Plot for quality of Data prediction





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

FILE_SAVE= paste0(FILE_OUTPUT_PLOTS,"/VlnPlot_prediction_score_AllNK_TangBlood" ,".png")

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p17)
dev.off()



#table(data.query$predicted.id, data.query$meta_tissue)
for (data_query_tissue in list_data.query_tissue){
  for(pop_predicted in levels(data_query_tissue$predicted.id)){
    data.query_subset = subset(data_query_tissue, subset = predicted.id == pop_predicted)
    
    data.query_subset@meta.data = droplevels(data.query_subset@meta.data)
    
    Vlnp_loop3= VlnPlot(data.query_subset , features= "prediction.score.max", group.by = "meta_histology_relabel" , pt.size = 0, sort = "decreasing", y.max = 1.1) & 
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





p17b = VlnPlot(data.query , features= "prediction.score.max", group.by = "meta_histology_relabel", pt.size = 0, sort = "decreasing", y.max = 1.1 ) & stat_summary(fun.data=data_summary,color="black") & ggtitle( "All NK") & plot_annotation(title = paste(unique(data.query$meta_tissue), collapse = " ") ) &
  theme(plot.title = element_text(hjust=0.5, size= 6))   &      theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45, size = 7), axis.text.y = element_text(hjust=1,  size = 7) )

cat("\n\n")
print(p17b)
cat("\n\n")


FILE_SAVE= paste0(FILE_OUTPUT_PLOTS,"/VlnPlot_prediction_score_AllNK_Blood" ,".png")

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p17b)
dev.off()





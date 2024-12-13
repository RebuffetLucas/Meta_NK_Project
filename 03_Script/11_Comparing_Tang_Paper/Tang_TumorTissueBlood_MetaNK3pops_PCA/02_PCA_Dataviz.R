#Only vizualisation of the PCA

#Load the table

#For tumor + blood:
mat= readRDS(paste0(PATH_ANALYSIS_OUTPUT, "/mat_all_tum_only.rds")) #All (Tumor + Blood)
TISSUE= "Tumor"

mat= readRDS(paste0(PATH_ANALYSIS_OUTPUT, "/mat_all_blood_Only.rds")) #All (Tumor + Blood)
TISSUE= "Blood"

mat= readRDS(paste0(PATH_ANALYSIS_OUTPUT, "/mat_all_tum_blood.rds")) #All (Tumor + Blood)
TISSUE= "Tumor_and_Blood"

#Run PCA with ade4

res.pca <- dudi.pca(mat,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5 ,          # Nombre d'axes gardés
                    center= TRUE ,
                    scale= TRUE
)


#Information carried by the PCAs

p_infoPCs = fviz_eig(res.pca, addlabels = TRUE, ylim = c(0,20 ))

print(p_infoPCs)

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/Information_pct_PCA.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_infoPCs)
dev.off()



#Composition of the PCs
  #Composition of the variables soleil all
p_soleil = fviz_pca_var(res.pca,
                        col.var = "cos2", # Color by contributions to the PC
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE,     # Avoid text overlapping,
                        #select.var = list(cos2 = 0.4 )
                        
)

sorted_cos2 <- sort(p_soleil$data$cos2, decreasing = TRUE)
cos2_Treshold =  sorted_cos2[NUMBER_PCA_COMPONENT_psoleil]


png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/Soleil_AllVar.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_soleil)
dev.off()



  #Composition of the variables soleil top

p_soleil2 = fviz_pca_var(res.pca,
                         col.var = "cos2", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE,     # Avoid text overlapping,
                         select.var = list(cos2 = cos2_Treshold)
                         
)

p_soleil2

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/Soleil_TOP", NUMBER_PCA_COMPONENT_psoleil ,".png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_soleil2)
dev.off()

  #Composition of the variables Barplots

for (compteur in 1:5){
  print(compteur)
  plot = fviz_cos2(res.pca, choice = "var", axes = compteur, top = NUMBER_VAR_BARPLOTS_PCA_COMPO )
  png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/VAR_CONTRIB_PCA","/PC",compteur, "_TOP", NUMBER_VAR_BARPLOTS_PCA_COMPO ,".png"), width = 30, height = 15,  units = "cm", res=600 )
  print(plot)
  dev.off()
}


#PCA visualization

plot_fviz = fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

data_frame_p3 = plot_fviz$data
split_list=str_split(rownames(data_frame_p3), "_", 3)
split_df = data.frame( do.call("rbind", split_list))

data_frame_p3$tissue = split_df$X1
data_frame_p3$patho = split_df$X2 
data_frame_p3$group = split_df$X3

data_frame_p3$clean_name= gsub(pattern = "Blood", replacement = "B", x = data_frame_p3$name  )
data_frame_p3$clean_name= gsub(pattern = "Tumor", replacement = "T", x = data_frame_p3$clean_name  )
data_frame_p3$clean_name= gsub(pattern = "_NK[0-9]", replacement = "", x = data_frame_p3$clean_name  )



#Make PCA plot colored by tissue
p3 = ggplot(data_frame_p3, aes(x = x, y= y, col = tissue, shape = group, label= clean_name))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  geom_text_repel(size = PCA_PLOT_TEXT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#DF6589FF", "#3C1053FF"))+ 
  labs(x = "PC 1", y = "PC 2")+
  ggtitle(TISSUE)
  
p3

p3bis = ggplot(data_frame_p3, aes(x = x, y= y, col = tissue, shape = group, label= NULL))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#DF6589FF", "#3C1053FF"))+ 
  labs(x = "PC 1", y = "PC 2")+
  ggtitle(TISSUE)

p3bis


png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/PCA_Plots", "/PCA_ColbyTissue.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p3)
dev.off()

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/PCA_Plots", "/PCA_ColbyTissue_NoName.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p3bis)
dev.off()

#Fig7b
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig7b.pdf",
  plot = p3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



#Make PCA plot colored by cluster
p4 = ggplot(data_frame_p3, aes(x = x, y= y, col = group, shape = tissue, label= clean_name))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  geom_text_repel(size = PCA_PLOT_TEXT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values =  c("#F8766D" , "#8494FF", "#0CB702" ))+ 
  labs(x = "PC 1", y = "PC 2")+
  ggtitle(TISSUE)

p4

p4bis = ggplot(data_frame_p3, aes(x = x, y= y, col = group, shape = tissue, label= NULL))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#F8766D" , "#8494FF", "#0CB702" ))+ 
  labs(x = "PC 1", y = "PC 2")+
  ggtitle(TISSUE)

p4bis


png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/PCA_Plots", "/PCA_ColbyCluster.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p4)
dev.off()

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/PCA_Plots", "/PCA_ColbyCluster_NoName.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p4bis)
dev.off()




#Covariance Heatmap


covMat = cor(t(mat), method = "spearman") #Spearman = nonparametric, pearson = parametric


#Fig7b
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig7b.pdf",
  plot = p3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 15,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


data_frame_heatmap = as.data.frame(covMat)
split_list_heatmap=str_split(rownames(data_frame_heatmap ), "_", 3)
split_df_heatmap = data.frame( do.call("rbind", split_list_heatmap))
colnames(split_df_heatmap) = c("Tissue", "Patho", "Cluster")


Pheatmap_plot  = pheatmap(covMat, border_color = "white", fontsize = 3, annotation_row = split_df_heatmap, annotation_colors = my_colourPheatmap)
Pheatmap_plot

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE , "/Heatmap_Covariance.png"), width = 20, height = 15,  units = "cm", res=600 )
print(Pheatmap_plot)
dev.off()



#Visualize PC3

#PCA visualization

plot_fviz_PC3 = fviz_pca_ind(res.pca, axes = c(2,3) ,col.ind = "cos2", 
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE # Avoid text overlapping (slow if many points)
)

data_frame_p3_PC3 = plot_fviz_PC3$data
split_list=str_split(rownames(data_frame_p3_PC3), "_", 3)
split_df = data.frame( do.call("rbind", split_list))

data_frame_p3_PC3$tissue = split_df$X1
data_frame_p3_PC3$patho = split_df$X2 
data_frame_p3_PC3$group = split_df$X3

data_frame_p3_PC3$clean_name= gsub(pattern = "Blood", replacement = "B", x = data_frame_p3_PC3$name  )
data_frame_p3_PC3$clean_name= gsub(pattern = "Tumor", replacement = "T", x = data_frame_p3_PC3$clean_name  )
data_frame_p3_PC3$clean_name= gsub(pattern = "_NK[0-9]", replacement = "", x = data_frame_p3_PC3$clean_name  )


p4_PC3 = ggplot(data_frame_p3_PC3, aes(x = x, y= y, col = group, shape = tissue, label= clean_name))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  geom_text_repel(size = PCA_PLOT_TEXT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values =  c("#F8766D" , "#8494FF", "#0CB702" ))+ 
  labs(x = "PC 2", y = "PC 3")+
  ggtitle(TISSUE)

p4_PC3

p4bis_PC3 = ggplot(data_frame_p3_PC3, aes(x = x, y= y, col = group, shape = tissue, label= NULL))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#F8766D" , "#8494FF", "#0CB702" ))+ 
  labs(x = "PC 2", y = "PC 3")+
  ggtitle(TISSUE)

p4bis_PC3


png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/PCA_Plots", "/PCA_ColbyCluster_PC3.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p4_PC3)
dev.off()

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/PCA_Plots", "/PCA_ColbyCluster_NoName_PC3.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p4bis_PC3)
dev.off()


#Fig7c
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig7c.pdf",
  plot = p4_PC3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 13,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#Fig7d
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig7d.pdf",
  plot = p4_PC3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 13,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



#Making  plots to compare PCA across clusters:

#Prepare the data frame
data_Frame_Violin = res.pca$li
data_Frame_Violin$name = rownames(data_Frame_Violin)

split_list=str_split(rownames(data_Frame_Violin), "_", 3)
split_df = data.frame( do.call("rbind", split_list))

data_Frame_Violin$tissue = split_df$X1
data_Frame_Violin$patho = split_df$X2 
data_Frame_Violin$group = split_df$X3

  #Drawing a ViolinPlot
for(i in 1:5) {
  axis_name <- paste("Axis", i, sep = "")
  p <- ggplot(data_Frame_Violin, aes_string(x = "group", y = axis_name, fill= "group")) + 
    geom_violin(trim = FALSE) + 
    labs(title = paste("Violin plot of", axis_name, "by group"), x = "Group", y = axis_name) +
    theme_bw()+
    scale_fill_manual(values = c("#F8766D" , "#8494FF", "#0CB702" ))
    
  print(p)
  
  png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_",TISSUE ,"/Violin_PCA_Values", "/Vln_PC", i, ".png"), width = 15, height = 10,  units = "cm", res=600 )
  print(p)
  dev.off()
  
}



  # Performing Kruskal-Wallis test for each Axis
for(i in 1:5) {
  axis_name <- paste("Axis", i, sep = "")
  print(paste("Kruskal-Wallis test for", axis_name))
  formula <- as.formula(paste(axis_name, "~ group"))
  test_result <- kruskal.test(formula, data = data_Frame_Violin)
  print(test_result)
}

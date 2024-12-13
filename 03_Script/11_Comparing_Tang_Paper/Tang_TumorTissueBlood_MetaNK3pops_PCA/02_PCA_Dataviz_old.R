#Only vizualisation of the PCA


#Load the table

mat= readRDS(paste0(PATH_ANALYSIS_OUTPUT, "/mat_all_tum_blood.rds")) #All (Tumor + Blood)


#Run PCA with ade4

res.pca <- dudi.pca(mat,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5 ,          # Nombre d'axes gardés
                    center= TRUE ,
                    scale= TRUE
)


#Basic
fviz_eig(res.pca)



df= as.data.frame( mat )
df$group = str_sub(rownames(df),-3,-1)
df$group = as.factor(df$group)

df$tissue = str_sub(rownames(df),1,5)
df$tissue = as.factor(df$tissue)


#PCA viz data colr
p_cluster = fviz_pca_ind(res.pca,
                         col.ind = df$group, 
                         palette = palette3[c(1,3,2)],
                         repel = TRUE ,
                         labelsize= 7,
                         pointsize= 2,
                         mean.point= FALSE
)  + theme(text = element_text(size = 12))  +   labs(title ="PCA", x = "PC1", y = "PC2")  



print(p_cluster)

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_No_Elypse.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_cluster)
dev.off()


p3 = ggplot(cbind(p_cluster$data,df[,c("group","tissue")]),
            aes(x=x,y=y,col=group,shape=tissue)) + 
  geom_point() + 
  theme_bw() +   
  labs(title ="PCA", x = "PC1", y = "PC2")  +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')

p_cluster$data$clean_name <- gsub(pattern = "_NK+[0-9]", replacement = "", x = p_cluster$data$name)
p_cluster$data$clean_name <- gsub(pattern = "Blood", replacement = "B", x = p_cluster$data$clean_name)
p_cluster$data$clean_name <- gsub(pattern = "Tumor", replacement = "T", x = p_cluster$data$clean_name)

df_plot_3 <- cbind(p_cluster$data,df[,c("group","tissue")])

p3 <- ggplot(df_plot_3, aes(x=x,y=y,col=group,shape=tissue, label =clean_name)) + 
  geom_point(size = 3) + 
  geom_text_repel(size = 4) +
  theme_bw() +   
  labs(title ="PCA", x = "PC1", y = "PC2")  +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')

p3_bis <- ggplot(df_plot_3, aes(x=x,y=y,col=group,shape=tissue)) + 
  geom_point(size = 2) + 
  theme_bw() +   
  labs(title ="PCA", x = "PC1", y = "PC2")  +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')

print(plot_grid(p3, p3_bis))


p_tissue = fviz_pca_ind(res.pca,
                        col.ind = df$tissue, 
                        palette = c("#F8766D" ,"#8494FF"),
                        repel = TRUE ,
                        labelsize= 7,
                        pointsize= 2,
                        mean.point= FALSE
)  + theme(text = element_text(size = 12))  +   labs(title ="PCA", x = "PC1", y = "PC2")  



print(p_tissue)

p4 = ggplot(cbind(p_tissue$data,df[,c("group","tissue")]),
            aes(x=x,y=y,col=group,shape=tissue)) + geom_point() + theme_classic() +   labs(title ="PCA", x = "PC1", y = "PC2")  



png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_No_Elypse.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_cluster)
dev.off()







#With Elypses

p1 = fviz_pca_ind(res.pca,
                  #col.ind = "coord", 
                  #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                  habillage=df$group,
                  axes= c(1,2),
                  cols.ind = palette3,
                  invisible='quali',
                  repel = TRUE ,
                  labelsize= 8,
                  pointsize= 3,
                  addEllipses= TRUE,
                  #ellipse.level=1.8,
                  legend= "none"
)  + theme(text = element_text(size = 20)) #+ ggforce::geom_mark_ellipse(aes(fill = df$group ,  color = df$group)) +  coord_equal()

print(p1)

png(file=paste0(FILE_OUTPUT_PLOTS , "/PCA_With_Elypse.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()












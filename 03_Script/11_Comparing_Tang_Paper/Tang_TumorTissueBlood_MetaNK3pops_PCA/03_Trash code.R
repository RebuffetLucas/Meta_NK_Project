# TRash code for data viz






#USELESS STUFF
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












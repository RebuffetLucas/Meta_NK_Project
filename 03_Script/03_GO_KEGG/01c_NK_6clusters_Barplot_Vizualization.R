#Simple vizualisation
#Not used in the pipeline

#Table of extraction
#For TableS3 page 1

lanes_of_interest= c(20,23,35,40,46,57,82,86) #Lanes of interest of unique(top10$Description)
Description_of_interest  = unique(top10$Description)[lanes_of_interest]








#  !!!!!!!    CODE à adapter pour eviter le soucis avec les descrptions redondantes !!!!!!!! #



dfm  = top10[top10$Description %in%Description_of_interest,]

dfm = top10

ggbarplot(dfm, x = "Description", y = "GeneRatio",
          fill = "Cluster",               # change fill color by cyl
          color = "Black",            # Set bar border colors to black
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal(),
          width= 0.7
) +  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + theme(text = element_text(size = 15))


ggbarplot(dfm, x = "Description", y = "GeneRatio",
          fill = "Cluster",               # change fill color by cyl
          color = "Black",            # Set bar border colors to black
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal(),
          width= 0.7
) +  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + theme(text = element_text(size = 15)) + coord_flip()


dfm$Cluster <- factor(dfm$Cluster, levels = c("NK_1", "NK_2", "NK_3"))

levels(dfm$Cluster)


#saveRDS(dfm,"/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Review_Nature/05_Output/DEG_NK1_2_3/best_genes_for_figures.rds")

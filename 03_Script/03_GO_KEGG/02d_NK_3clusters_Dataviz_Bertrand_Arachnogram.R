#Bertrand Fancy vizualization

lanes_of_interest= c(23,43,50,55,66,69,106,107) #Lanes of interest of unique(top10$Description)

GO_to_plot  = unique(top10$Description)[lanes_of_interest]

GO_to_plot = c("regulation of cytoskeleton organization"    ,     "immune response-activating cell surface receptor signaling pathway",
"positive regulation of leukocyte activation"    ,                "cell killing"   ,                                                   
"antimicrobial humoral response"     ,   "positive regulation of leukocyte cell-cell adhesion"   ,            
"regulation of leukocyte chemotaxis"     ,                            "regulation of leukocyte differentiation"        )                   
 

Matrix_Compare <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enrichGO", OrgDb = "org.Hs.eg.db",ont= "BP", universe = UNIVERSE_EntrezID, readable = TRUE, pvalueCutoff = 1, qvalue=1)



Matrix_Compare = setReadable(Matrix_Compare , OrgDb = org.Hs.eg.db, keyType="ENTREZID")
Matrix_Compare2 = Matrix_Compare@compareClusterResult
Matrix_Compare3  = Matrix_Compare2[Matrix_Compare2$Description %in% GO_to_plot, ]


MT = data.frame(Matrix_Compare3)
MT$GeneRatio= parse_ratio(MT$GeneRatio)
MT2= MT[,c("Cluster", "Description", "p.adjust")]

MT2$GeneRatio = -log10(MT2$p.adjust)

MT2$GeneRatio
table(MT2$Description)

# Make the plot
p = ggplot(MT2, aes(x=as.factor(Description), y=GeneRatio, fill=Cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=1, position="dodge") +
  
  
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = 0.45, y = c(2, 4,6), label = c("2", "4", "6") , color="grey20", size=3 , angle=0, fontface="bold", hjust=1) +
  annotate("text", x = 0.45, y = 6.5, label = "-log10(p.adj)" , color="grey20", size=3 , angle=0, fontface="bold", hjust=1) +
  
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=0.5, position="dodge") +
  ylim(-0.25,7)+
  theme_minimal() +
  theme(
    axis.text.x = element_text(size=6,face="bold",vjust=-100),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),legend.position="bottom",legend.justification = c(1,1),legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid",colour ="black"),
    ) + geom_hline(yintercept = -log10(0.05),linetype="dashed")+geom_hline(yintercept = c(2,4,6),linetype="dashed",color = "grey20",alpha=0.3)+
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5),linetype="dashed",color = "black",alpha=0.3)+
  coord_polar(start=0)
# Add base line information
# geom_segment(aes(x = 0, y = -0.5, xend = 19, yend = -0.5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )

p

# Make the plot without legend
p2 = ggplot(MT2, aes(x=as.factor(Description), y=GeneRatio, fill=Cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=1, position="dodge") +
  
  
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = 0.45, y = c(2, 4,6), label = c("2", "4", "6") , color="grey20", size=4 , angle=0, fontface="bold", hjust=1) +
  #annotate("text", x = 0.45, y = 8.5, label = "-log10(p.adj)" , color="grey20", size=4 , angle=0, fontface="bold", hjust=1) +
  
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=0.5, position="dodge") +
  ylim(-0.25,7)+
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),legend.position="bottom",legend.justification = c(1,1),legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid",colour ="black"),
  ) + geom_hline(yintercept = -log10(0.05),linetype="dashed")+geom_hline(yintercept = c(2,4,6),linetype="dashed",color = "grey20",alpha=0.3)+
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5),linetype="dashed",color = "black",alpha=0.3)+
  coord_polar(start=0)
# Add base line information
# geom_segment(aes(x = 0, y = -0.5, xend = 19, yend = -0.5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )

p2



png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures_V3/3clusters_WithLegend.png", width = 20, height =20, units= "cm" ,  res=600 )
p
dev.off()


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures_V3/3clusters_NoLegend.png", width = 20, height =20, units= "cm" ,  res=600 )
p2
dev.off()






pdf("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/3clusters_NoLegend.pdf")
p2
dev.off()

svg("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/3clusters_NoLegend.svg")
p2
dev.off()



pdf("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/3clusters_WithLegend.pdf")
p
dev.off()

svg("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/3clusters_WithLegend.svg")
p
dev.off()



#Bertrand Fancy vizualization on Reactome shit

#Reactome
"Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins."
"The citric acid (TCA) cycle and respiratory electron transport" 
"Respiratory electron transport"    
"Negative regulation of NOTCH4 signaling" 
"Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha"
"Regulation of Apoptosis"                                                                                            
"Cellular response to hypoxia"
"Metabolism of polyamines"  
"Interleukin-12 signaling"  
"Regulation of mitotic cell cycle"
"Signaling by NOTCH"
"FCGR3A-mediated phagocytosis" 



GO_to_plot = c(
  "Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.",
  "The citric acid (TCA) cycle and respiratory electron transport" ,
  "Respiratory electron transport"    ,
  "Negative regulation of NOTCH4 signaling", 
  "Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha",
  "Regulation of Apoptosis"  ,                                                                                          
  "Cellular response to hypoxia",
  "Metabolism of polyamines"  ,
  "Interleukin-12 signaling"  ,
  "Regulation of mitotic cell cycle",
  "Signaling by NOTCH",
  "FCGR3A-mediated phagocytosis" 
)

Matrix_Compare <- compareCluster(ENTREZID~cluster, data=ORA_DF, fun="enrichPathway", universe = UNIVERSE_EntrezID, readable = TRUE, pvalueCutoff = 1, qvalue=1)


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
  annotate("text", x = 0.45, y = c(2, 4,6,8), label = c("2", "4", "6", "8") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  annotate("text", x = 0.45, y = 8.5, label = "-log10(FDR)" , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=0.5, position="dodge") +
  ylim(-0.25,9)+
  theme_minimal() +
  theme(
    axis.text.x = element_text(size=6,face="bold",vjust=-100),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),legend.position="bottom",legend.justification = c(1,1),legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid",colour ="black"),
    ) + geom_hline(yintercept = -log10(0.05),linetype="dashed")+geom_hline(yintercept = c(2,4,6,8),linetype="dashed",color = "grey",alpha=0.3)+
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5),linetype="dashed",color = "grey",alpha=0.3)+
  coord_polar(start=0)
# Add base line information
# geom_segment(aes(x = 0, y = -0.5, xend = 19, yend = -0.5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )

p

# Make the plot without legend
p2 = ggplot(MT2, aes(x=as.factor(Description), y=GeneRatio, fill=Cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=1, position="dodge") +
  
  
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = 0.45, y = c(2, 4,6,8), label = c("2", "4", "6", "8") , color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
  #annotate("text", x = 0.45, y = 8.5, label = "-log10(FDR)" , color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
  
  
  geom_bar(aes(x=as.factor(Description), y=GeneRatio, fill=Cluster), stat="identity", alpha=0.5, position="dodge") +
  ylim(-0.25,9)+
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),legend.position="bottom",legend.justification = c(1,1),legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid",colour ="black"),
  ) + geom_hline(yintercept = -log10(0.05),linetype="dashed")+geom_hline(yintercept = c(2,4,6,8),linetype="dashed",color = "grey",alpha=0.3)+
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5),linetype="dashed",color = "grey",alpha=0.3)+
  coord_polar(start=0)
# Add base line information
# geom_segment(aes(x = 0, y = -0.5, xend = 19, yend = -0.5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )

p2


pdf("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/6clusters_NoLegend.pdf")
p2
dev.off()

svg("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/6clusters_NoLegend.svg")
p2
dev.off()



pdf("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/6clusters_WithLegend.pdf")
p
dev.off()

svg("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/03_GO_KEGG/Figures/6clusters_WithLegend.svg")
p
dev.off()



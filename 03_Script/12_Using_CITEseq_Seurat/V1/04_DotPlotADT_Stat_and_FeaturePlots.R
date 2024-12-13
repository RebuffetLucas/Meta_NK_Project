#Try sub-clustering of the CITE-seq Dataset
PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Review_Nature/05_Output/CITEseq_Analysis/Seurat_3Clusters.rds")


REMOVE = c("CD38-2", "CD44-2", "CD56-2", "CD158") #Remove the ADT that are present twice or undefined
PBMC  = subset(PBMC , subset = time == "0")
table(PBMC$wsnn_res.0.2 )
levels(PBMC$wsnn_res.0.2)  = c("NK1", "NK3", "NK2")

PBMC$wsnn_res.0.2  = factor(PBMC$wsnn_res.0.2, levels = c("NK1", "NK2", "NK3"))

PBMC$wsnn_res.0.2 <-factor(PBMC$wsnn_res.0.2 ,c("NK1", "NK2", "NK3"))

PBMC= SetIdent(PBMC, value="wsnn_res.0.2" )

DimPlot(PBMC, reduction = "wnn.umap")


#Get the list of genes present in at list X percent of the cells
MIN_PCT_PROT= 20
NUMBER_TOP_SCORING = 20
palette<-c('NK1'='#F8766D','NK2'='#8494FF',
           'NK3'='#0CB702')

DimPlot(PBMC, reduction = "wnn.umap", cols = palette)



#Normalize both ADT and RNA
#Normalize RNA expression
DefaultAssay(PBMC) <- 'SCT'
PBMC= SCTransform(PBMC, assay= "SCT", )

#Normalize Protein expression
DefaultAssay(PBMC) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(PBMC) <- rownames(PBMC[["ADT"]])
PBMC <- NormalizeData(PBMC, normalization.method = 'CLR', margin = 2)


PBMC= SetIdent(PBMC, value = "wsnn_res.0.2")
PBMC$wsnn_res.0.2= droplevels(PBMC$wsnn_res.0.2)



#Have a look at markers
#PBMC = subset(Seurat_NK, subset= time=="0" )

p0 = DimPlot(PBMC, reduction = "wnn.umap", cols = palette) +NoAxes()

p0


Markers_Prots = FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

markers_overCl=Markers_Prots

topGenes=c()
for (cl in sort(unique(Idents(PBMC)))) {
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
  print(	topGenes_cl)
}

topGenes2=topGenes

topGenesb=unique(unlist(topGenes))

topGenesb = topGenesb[ ! topGenesb %in% REMOVE]



#Visualisation
DotPlot(PBMC, features = topGenesb , cols = "Spectral", scale=FALSE)  + coord_flip()
p_dotplot2 = DotPlot(PBMC, features = topGenesb , cols = "RdBu")  + coord_flip()
p_dotplot = DotPlot(PBMC, features = topGenesb , cols = "Spectral")  + coord_flip()


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3clusters_DotPlot_ADT.png" , width = 12, height = 20,  units = "cm", res=600 )
p_dotplot
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3clusters_DotPlot_ADT_RdBu.png" , width = 12, height = 20,  units = "cm", res=600 )
p_dotplot2
dev.off()

#Save as pdf figures for NI


pdf(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3clusters_DotPlot_ADT_RdBu.pdf" , paper ="a4" , width = 4.5, height = 6.8 )
p_dotplot2
dev.off()

#Fig1a
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig1a.pdf",
  plot = p0,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 18,
  height = 15,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#Fig1c
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig1c.pdf",
  plot = p_dotplot2,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 12,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)








#Feature Plots



p0
p1 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[1:4], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p2 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[5:8], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p3 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[9:12], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p4 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[13:16] , min.cutoff = 'q01', max.cutoff = 'q99') &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p5 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[17:20] , min.cutoff = 'q01', max.cutoff = 'q99') &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p6 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[21:24], min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p7 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD56-1", "CD16", "CD38-1", "CD57") , min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

p8 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD56-1", "CD16", "CD38-1", "CD122") , min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p9 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD57", "CD44-1", "CD27", "CX3CR1"), min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p10 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD54", "CD45RB", "CD335", "CD96"), min.cutoff = 'q01', max.cutoff = 'q99' ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()


p10

summary(PBMC@assays[["ADT"]]@data['CD57',])



#p1= cowplot::plot_grid( p0, p1a, p1b, p1c, ncol=2)

p2
p3
p4
p5
p6
p7



png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P1.png" , width = 20, height = 20,  units = "cm", res=600 )
p1
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P2.png" , width = 20, height = 20,  units = "cm", res=600 )
p2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P3.png" , width = 20, height = 20,  units = "cm", res=600 )
p3
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P4.png" , width = 20, height = 20,  units = "cm", res=600 )
p4
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P5.png" , width = 20, height = 20,  units = "cm", res=600 )
p5
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P6.png" , width = 20, height = 20,  units = "cm", res=600 )
p6
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P7_56_16_38_57.png" , width = 20, height = 20,  units = "cm", res=600 )
p7
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P8_16_56_38_122.png" , width = 20, height = 20,  units = "cm", res=600 )
p8
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P9_57_44_27_CX3CR1.png" , width = 20, height = 20,  units = "cm", res=600 )
p9
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_P10_54_96_335_45.png" , width = 20, height = 20,  units = "cm", res=600 )
p10
dev.off()


#Version gradients de rouge
p1 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[1:4], min.cutoff = 'q01', max.cutoff = 'q99', pt.size = 0.8 ) &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p2 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[5:8], min.cutoff = 'q01', max.cutoff = 'q99', pt.size = 0.8 ) &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p3 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[9:12], min.cutoff = 'q01', max.cutoff = 'q99', pt.size = 0.8 ) &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p4 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[13:16] , min.cutoff = 'q01', max.cutoff = 'q99', pt.size = 0.8)  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p5 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[17:20] , min.cutoff = 'q01', max.cutoff = 'q99', pt.size = 0.8)  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p6 = FeaturePlot(PBMC, reduction= "wnn.umap", features= topGenesb[21:24], min.cutoff = 'q01', max.cutoff = 'q99' , pt.size = 0.8)  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p7 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD56-1", "CD16", "CD38-1", "CD57") , min.cutoff = 'q01', max.cutoff = 'q99', pt.size = 0.8 )  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()


p8 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD56-1", "CD16", "CD38-1", "CD122") , min.cutoff = 'q01', max.cutoff = 'q99' , pt.size = 0.8)  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p9 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD57", "CD44-1", "CD27", "CX3CR1"), min.cutoff = 'q01', max.cutoff = 'q99' , pt.size = 0.8)  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
p10 = FeaturePlot(PBMC, reduction= "wnn.umap", features= c("CD54", "CD45RB", "CD335", "CD96"), min.cutoff = 'q01', max.cutoff = 'q99' , pt.size = 0.8)  &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()



p1
p2
p3
p4
p5
p6
p7
p8
p9


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/WNN_UMAP.png" , width = 20, height = 20,  units = "cm", res=600 )
p0
dev.off()


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP1.png" , width = 20, height = 20,  units = "cm", res=600 )
p1
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP2.png" , width = 20, height = 20,  units = "cm", res=600 )
p2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP3.png" , width = 20, height = 20,  units = "cm", res=600 )
p3
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP4.png" , width = 20, height = 20,  units = "cm", res=600 )
p4
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP5.png" , width = 20, height = 20,  units = "cm", res=600 )
p5
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP6.png" , width = 20, height = 20,  units = "cm", res=600 )
p6
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP7_56_16_38_57.png" , width = 20, height = 20,  units = "cm", res=600 )
p7
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP8_16_56_38_122.png" , width = 20, height = 20,  units = "cm", res=600 )
p8
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP9_57_44_27_CX3CR1.png" , width = 20, height = 20,  units = "cm", res=600 )
p9
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/FeaturePlot_ADT_RedsVersionP10_54_96_335_45.png" , width = 20, height = 20,  units = "cm", res=600 )
p10
dev.off()


#Fig1d1
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig1d1.pdf",
  plot = p8,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#Fig1d2
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig1d2.pdf",
  plot = p9,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#Look at signatures
#have a look at their signatures
DefaultAssay(object = PBMC) <- "SCT"
PBMC= SCTransform(PBMC, assay= "SCT")



All_Markers = FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All


markers_overCl=top_All

topGenes=c()
for (cl in sort(unique(Idents(PBMC)))) {
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
  print(	topGenes_cl)
}

topGenes2=topGenes

topGenesb=unique(unlist(topGenes))

p_dotplot3 = DotPlot(PBMC, features = topGenesb , cols = "Spectral")  + coord_flip()
p_dotplot4 = DotPlot(PBMC, features = topGenesb , cols = "RdBu")  + coord_flip()


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3clusters_DotPlot_RNA.png" , width = 12, height = 35,  units = "cm", res=600 )
p_dotplot3
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3clusters_DotPlot_RNA_RdBu.png" , width = 12, height = 35,  units = "cm", res=600 )
p_dotplot4
dev.off()


#Fig1b
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig1b.pdf",
  plot = p_dotplot4,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 35,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)





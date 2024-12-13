#Try sub-clustering of the CITE-seq Dataset
PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Review_Nature/05_Output/CITEseq_Analysis/Seurat_3Clusters.rds")

PBMC  = subset(PBMC , subset = time == "0")
table(PBMC$wsnn_res.0.2 )
levels(PBMC$wsnn_res.0.2 )
PBMC$wsnn_res.0.2 = relevel(PBMC$wsnn_res.0.2, "NK_2")
PBMC$wsnn_res.0.2 = relevel(PBMC$wsnn_res.0.2, "NK_1")

DimPlot(PBMC, reduction = "wnn.umap")
PBMC= SetIdent(PBMC, value="wsnn_res.0.2" )
DimPlot(PBMC, reduction = "wnn.umap", group.by= "celltype.l3") 

#Have a look at markers
#PBMC = subset(Seurat_NK, subset= time=="0" )

p1 = DimPlot(PBMC, reduction = "wnn.umap") +NoAxes()

p1


#Get the list of genes present in at list X percent of the cells
MIN_PCT_PROT= 20

#Normalize both ADT and RNA
  #Normalize RNA expression
DefaultAssay(PBMC) <- 'SCT'
PBMC= SCTransform(PBMC, assay= "SCT")

  #Normalize Protein expression
DefaultAssay(PBMC) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(PBMC) <- rownames(PBMC[["ADT"]])
PBMC <- NormalizeData(PBMC, normalization.method = 'CLR', margin = 2)



PBMC= SetIdent(PBMC, value = "wsnn_res.0.2")
PBMC$wsnn_res.0.2= droplevels(PBMC$wsnn_res.0.2)


#HeatMap RNA level
DefaultAssay(PBMC) <- 'ADT'

#Heatmap List
List_Genes1 = c("CD45RA", "CD57", "CD16",  "CD2", "CD56-1",  "CD161", "CD38-1") #Maturation
List_Genes2 = c("CD319" ,  "CD337","CD314" , "CD59", "CD137" , "CD335", "CD226") #Activating Receptors
List_Genes3 = c("CD122",  "CD127", "CD25", "CD195") #IL receptors
List_Genes4 = c( "CD366","CD158e1","CD158", "CD158b",  "CD158f" , "CD223", "TIGIT", "CD96", "CD279") # Récepteurs inhibiteurs
List_Genes5 = c(  "CD107a" ,"CD244","CD69","CD103") #Other markers of interest ADD    CD49a
List_Genes6 = c(  "CD56-2","CD195", "CD25", "CD158f" , "CD279" , "CD49a",   "CD226") #Not in RNA

#New List

List_Genes1 = c("CD152", "CD223", "CD279", "TIGIT", "CD366") #Exhaustion
List_Genes2 = c("CD314", "CD335", "CD337" ,"CD54", "CD184", "CD117", "CD271", "CD178") #NK receptors
List_Genes3 = c("CD3-1", "CD3-2") #T like receptors
List_Genes4 = c("CD69", "CD103", "CD186") # Tissue Homing
List_Genes5 = c("CD158", "CD158b", "CD158e1", "CD158f") # KIR family
List_Genes6 = c("CD49a", "CD49b", "CD49d") #CD49 family


#intersect(All_genes , c("B3GAT1", "CD244", "CD25", "IL2RB", "IL2RG", "SIGLEC7"))

List_Markers = list.append(List_Genes1, List_Genes2, List_Genes3, List_Genes4 , List_Genes5  , List_Genes6)
List_of_Lists= list(List_Genes1, List_Genes2, List_Genes3, List_Genes4 , List_Genes5,  List_Genes6)

List_total= intersect(List_Markers, rownames(PBMC@assays[["ADT"]]))

length(List_total)

#Keep Genes expressed in more than X percent cells
perc_exp <- DotPlot(PBMC, features=List_total, group.by="wsnn_res.0.2")$data[, c("features.plot", "id", "pct.exp")]

perc_exp %>%
  filter(pct.exp > MIN_PCT_PROT ) -> Acceptable_Genes

Acceptable_Genes = unique(Acceptable_Genes$features.plot)
length(Acceptable_Genes)

setdiff(List_total, Acceptable_Genes)

List_total= intersect(List_Markers, Acceptable_Genes)
List_of_Lists=  lapply(List_of_Lists , function(x) intersect(x, Acceptable_Genes))


FeaturePlot(PBMC, features = List_Genes1, reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = List_Genes2, reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = List_Genes3, reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = List_Genes4, reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = List_Genes5, reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = List_Genes6, reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

VlnPlot(PBMC, features= List_Genes1, pt.size = 0 )
VlnPlot(PBMC, features= List_Genes2, pt.size = 0 )
VlnPlot(PBMC, features= List_Genes3, pt.size = 0 )
VlnPlot(PBMC, features= List_Genes4, pt.size = 0 )
VlnPlot(PBMC, features= List_Genes5, pt.size = 0 )
VlnPlot(PBMC, features= List_Genes6, pt.size = 0 )


p7 = FeaturePlot(PBMC, features = c("CD57"), reduction= "wnn.umap", max.cutoff = 4 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p8 = FeaturePlot(PBMC, features = c("CD56-1"), reduction= "wnn.umap", max.cutoff = 3 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p9 = FeaturePlot(PBMC, features = c("CD16"), reduction= "wnn.umap" ,  max.cutoff = 5)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

p1 + p7 + p8 + p9

ggsave("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3Clusters_wnnUMAPandFeaturePlotProt.png",
       plot = p1 + p7 + p8 + p9 , width = 20, height = 15, dpi = 600, units = "cm")



p10 = FeaturePlot(PBMC, features = c("CD161"), reduction= "wnn.umap", max.cutoff = 1.5 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p11 = FeaturePlot(PBMC, features = c("CD56-2"), reduction= "wnn.umap", max.cutoff = 4 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p12 = FeaturePlot(PBMC, features = c("CD38-1"), reduction= "wnn.umap", max.cutoff = 4)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

p10 + p11 + p12

p13 = FeaturePlot(PBMC, features = c("CD38-2"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p14 = FeaturePlot(PBMC, features = c("CD52"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()
p15 = FeaturePlot(PBMC, features = c("CD314"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()


p13 + p14 + p15

p16 = FeaturePlot(PBMC, features = c("CD314"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p17 = FeaturePlot(PBMC, features = c("CD44-1"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p18 = FeaturePlot(PBMC, features = c("CD54"), reduction= "wnn.umap" ,  max.cutoff = 2.1)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p16 + p17 + p18


p19 = FeaturePlot(PBMC, features = c("CD44-2"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p20 = FeaturePlot(PBMC, features = c("CD158"), reduction= "wnn.umap",  max.cutoff = 2)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p21 = FeaturePlot(PBMC, features = c("CD158b"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p19 + p20 + p21

p22 = FeaturePlot(PBMC, features = c("CD158e1"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))







p7 = FeaturePlot(PBMC, features = c("CD335"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD184"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD178"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD337"), reduction= "wnn.umap", max.cutoff = 2.5)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD117"), reduction= "wnn.umap", max.cutoff = 3)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD178"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD69"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD103"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD186"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD107a"), reduction= "wnn.umap", max.cutoff = 2.5)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD107a"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD107a"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p2 = FeaturePlot(PBMC, features = c("CD314"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p7 = FeaturePlot(PBMC, features = c("CD158"), reduction= "wnn.umap",  max.cutoff = 2)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD158b"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD158e1"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p2 + p7
p8 + p9

ggsave("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3Clusters_wnnUMAPandFeaturePlot_CD314_CD158.png",
       plot = p2 + p7 , width = 20, height = 15, dpi = 600, units = "cm")

ggsave("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/Final_Figures/CITEseq_3Clusters_wnnUMAPandFeaturePlot_CD158b_CD158e1.png",
       plot = p8 + p9 , width = 20, height = 15, dpi = 600, units = "cm")


p7 = FeaturePlot(PBMC, features = c("CD158f"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD158f"), reduction= "wnn.umap",  max.cutoff = 1)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD158f"), reduction= "wnn.umap",  max.cutoff = 0.5) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD158f"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD158f"), reduction= "wnn.umap",  max.cutoff = 1)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD158f"), reduction= "wnn.umap",  max.cutoff = 0.5) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD49a"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD49b"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD49d"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD152"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD223"), reduction= "wnn.umap", max.cutoff = 1.2)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD279"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9


p7 = FeaturePlot(PBMC, features = c("TIGIT"), reduction= "wnn.umap", max.cutoff = 1.5 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD366"), reduction= "wnn.umap", max.cutoff = 1.9)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD96"), reduction= "wnn.umap", max.cutoff = 1.4) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

List_Genes2 = c("CD319" , "CD59", "CD137" , "CD226")

p7 = FeaturePlot(PBMC, features = c("CD319"), reduction= "wnn.umap" )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD59"), reduction= "wnn.umap" )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD137"), reduction= "wnn.umap" ) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9


p7 = FeaturePlot(PBMC, features = c("CD226"), reduction= "wnn.umap", max.cutoff = 1.5 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD226"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD226"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9


#c("CD122",  "CD127", "CD25", "CD195") #IL receptors

p7 = FeaturePlot(PBMC, features = c("CD122"), reduction= "wnn.umap" )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD127"), reduction= "wnn.umap")  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD25"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9

p7 = FeaturePlot(PBMC, features = c("CD195"), reduction= "wnn.umap" )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p8 = FeaturePlot(PBMC, features = c("CD195"), reduction= "wnn.umap", max.cutoff = 2)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p9 = FeaturePlot(PBMC, features = c("CD195"), reduction= "wnn.umap") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p7 + p8 + p9





p10

FeaturePlot(PBMC, features = c("CD178"), reduction= "wnn.umap", max.cutoff = 2 )  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


FeaturePlot(PBMC, features = List_Genes2, reduction= "wnn.umap")


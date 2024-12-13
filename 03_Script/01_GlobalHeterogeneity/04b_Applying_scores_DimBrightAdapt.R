#Apply scores of Dim, Bright and Adapt to the initial UMAP

#All

data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}



#PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2_VF1_AllGenes_NewNames.rds")

PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

PBMC= SetIdent(PBMC, value= "FirstClust")
DimPlot(PBMC, label = TRUE, cols = palette)

PBMC= SetIdent(PBMC, value= "FirstClust")

PBMC$FirstClust = factor(PBMC$FirstClust, levels= ORDER_CLUST_LEGEND)

levels(PBMC$FirstClust)
levels(PBMC$SecondClust) = c("NK1C" ,  "NK3A" ,"NK1A", "NK3C", "NK1B", "NKint", "NK2" ,"NK3B")
PBMC$SecondClust = factor(PBMC$SecondClust, levels= ORDER_CLUST_LEGEND2)

DimPlot(PBMC, group.by = "SecondClust",label = TRUE, cols = palette2)


saveRDS(PBMC, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds" )


PBMC$SecondClust = droplevels(PBMC$SecondClust)


PBMC= subset(PBMC, idents= c("NK2A", "NK2B"))
PBMC= subset(PBMC, idents= c("NK2A", "NK2B"), invert= TRUE)

PBMC$FirstClust= droplevels(PBMC$FirstClust)

PBMC$SecondClust= droplevels(PBMC$FirstClust)


#Import signatures
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Crinier_13_NKgenes.xlsx") #13 genes Crinier 2018
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Ref_Crinier_2020.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Davis_Scores.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Davis_Scores_titledcluster.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Davis_Scores_titledclusterV2.xlsx")



FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Functionnal_Scores.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Zhang_Dim.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Zhang_Bright.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Zhang_BrightandDim.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Bhandoola.xlsx")

FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Zhang_BrightV2.xlsx")
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Zhang_DimV2.xlsx")

FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Zhang_DimV2.xlsx")


###Markers for scoring of subsets
# getting data from sheets
sheets_names <- openxlsx::getSheetNames(FILE_SIGNATURES_PATH)



#def function

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}



#For Malmberg
MONITORED_Markers = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg.rds")

MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)

MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(PBMC)))
MONITORED_Markers= lapply(MONITORED_Markers, FUN=function(x) head(x, n=40))

#Add the module scores
for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}





##### Vizualization #####

#Pheatmap

table = PBMC@meta.data[,c( "FirstClust" , paste0(names(MONITORED_Markers),1) ) ]

FeaturePlot(PBMC, features= "Circulatory_Score_GG1", pt.size = 0.8)

p1 = FeaturePlot(PBMC, features= "Skin_signature_GG1", pt.size = 0.6) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
p2= FeaturePlot(PBMC, features= "Circulatory_Score_GG1", pt.size = 0.6) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) &NoAxes()
p3 = FeaturePlot(PBMC, features= "Tissue_Residency_Schuster1", pt.size = 0.6) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
p4 =FeaturePlot(PBMC, features= "Tissue_Residency_Schuster_withF1", pt.size = 0.6) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()

p5 = FeaturePlot(PBMC, features= c("Skin_signature_GG1" , "Circulatory_Score_GG1" , "Tissue_Residency_Schuster1" , "Tissue_Residency_Schuster_withF1" ), pt.size = 0.6) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_Functionnal/FeaturePlot_Skin_signature_GG1.png", width = 18, height = 18,  units = "cm", res=600 )
p1
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_Functionnal/FeaturePlot_Circulatory_Score_GG1.png", width = 18, height = 18,  units = "cm", res=600 )
p2
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_Functionnal/FeaturePlot_Tissue_Residency_Schuster1.png", width = 18, height = 18,  units = "cm", res=600 )
p3
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_Functionnal/FeaturePlot_Tissue_Residency_Schuster_withF1.png", width = 18, height = 18,  units = "cm", res=600 )
p4
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_Functionnal/FeaturePlot_All_Four_Scoring.png", width = 22, height = 18,  units = "cm", res=600 )
p5
dev.off()



VlnPlot(PBMC, features= "Crinier_13_hNKGenes_20181" , pt.size = 0) & ggtitle("Score 13 NK Genes") & stat_summary(fun.data=data_summary,color="black")


VlnPlot(PBMC, features= "Tissue_Residency_Schuster1"  , pt.size = 0 ) & ggtitle("Tissue_Residency_Schuster1s") & stat_summary(fun.data=data_summary,color="black")






p3 = FeaturePlot(PBMC , features = paste0(names(MONITORED_Markers),1) , pt.size = 0.8) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) )

p1= VlnPlot(PBMC , features = paste0(names(MONITORED_Markers),1), group.by= "SecondClust", pt.size = 0, cols= palette2 ) 
p2 = VlnPlot(PBMC , features = paste0(names(MONITORED_Markers),1), group.by= "FirstClust" , pt.size = 0, cols = palette) & ggtitle("Score 13 NK Genes") & stat_summary(fun.data=data_summary,color="black")



png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_13Genes_NK/VlnPlot_13NKgenes_FirstClust.png", width = 18, height = 18,  units = "cm", res=600 )
p1
dev.off()

png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Scoring_13Genes_NK/VlnPlot_13NKgenes_SecondClust.png", width = 18, height = 18,  units = "cm", res=600 )
p2
dev.off()






df = table %>% 
  group_by(FirstClust) %>% 
  dplyr::summarise(across(everything(), list(mean)))

df= as.data.frame(df)
rownames(df)=df$FirstClust

df2= df[,-1]


pheatmap(as.matrix(df2), scale= "column", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_col = 15, fontsize_row = 15)

pheatmap(as.matrix(df2), scale= "column")


#VlnPlot

for (i in names(MONITORED_Markers)){
  print(VlnPlot(PBMC, features= paste0( i , "1") , pt.size = 0))
}

p1 = VlnPlot(PBMC, features = "NK_01" , pt.size = 0)
p2 = VlnPlot(PBMC, features = "NK_11" , pt.size = 0)
p3 = VlnPlot(PBMC, features = "NK_21", pt.size = 0)
p4 = VlnPlot(PBMC, features = "Adaptative_NK1", pt.size = 0)

p1 + p2+ p3 +p4
(p1 + p2 )
(p3 + p4 )

FeaturePlot(PBMC, feature = "ENKP_NK1" , pt.size = 0.3) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.4, 1.2) )
FeaturePlot(PBMC, feature = "ILCP_NK1" , pt.size = 0.3 )  +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.2, 0.5) )

FeaturePlot(PBMC, feature = "ENKP_NK1" , pt.size = 0.3 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) )
FeaturePlot(PBMC, feature = "ILCP_NK1" , pt.size = 0.3 )  +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) )

VlnPlot(PBMC, feature = "ENKP_NK1" , pt.size = 0.3 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) )
VlnPlot(PBMC, feature = "ILCP_NK1" , pt.size = 0.3 )  +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) )



+    scale_colour_gradient2( midpoint = median(PBMC$hNK_Bl11 )   , low= "blue", high = "red" , mid= "white" )

FeaturePlot(PBMC, feature = "hNK_Bl21" )  

FeaturePlot(PBMC, feature = c("GZMH", "CCL5", "KLRC2", "IL32") )  








#Feature Plot

#Inflam
FeaturePlot(PBMC, features = "Zhang_Inflammatory1", pt.size = 0.8 )  +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = "Zhang_Inflammatory1", pt.size = 0.8 , cols = viridis_plasma_dark_high) 
FeaturePlot(PBMC, features = "Zhang_Inflammatory1",  pt.size = 0.8 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-1, 1.8) )


FeaturePlot(PBMC, features = "Zhang_Inflammatory1" , pt.size = 1 ) +     scale_colour_gradient2( midpoint = median(PBMC$Zhang_Inflammatory1)   , low= "blue", high = "red" , mid= "white" )

#Cytotox
FeaturePlot(PBMC, features = "Zhang_Cytotoxicity1", pt.size = 0.6 )  +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(PBMC, features = "Zhang_Cytotoxicity1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-1, 2.4) )
FeaturePlot(PBMC, features = "Zhang_Cytotoxicity1" , pt.size = 0.6 ) +     scale_colour_gradient2( midpoint = median(PBMC$Zhang_Inflammatory1)   , low= "blue", high = "red" , mid= "white" )

#Circulatory_GG for Bright
FeaturePlot(PBMC, features = "Circulatory_Score_GG1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(0, 0.35) )
FeaturePlot(PBMC, features = "Skin_signature_GG1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(0, 0.6) )
FeaturePlot(PBMC, features = "Tissue_Residency_Schuster1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.5, 1) )
FeaturePlot(PBMC, features = "Tissue_Residency_Schuster_21",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.3, 0.5) )
FeaturePlot(PBMC, features = "TRM_Specific_Signature_Chen1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.5, 1) )

#Circulatory_GG for Dim
FeaturePlot(PBMC, features = "Circulatory_Score_GG1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.1, 0.32) )
FeaturePlot(PBMC, features = "Skin_signature_GG1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(0, 0.6) )
FeaturePlot(PBMC, features = "Tissue_Residency_Schuster1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.5, 1.1) )
FeaturePlot(PBMC, features = "Tissue_Residency_Schuster_21",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.3, 0.5) )
FeaturePlot(PBMC, features = "TRM_Specific_Signature_Chen1",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(-0.5, 0.9) )


FeaturePlot(PBMC, features = "GZMB",  pt.size = 0.6 ) +     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), limits= c(0, 3.5) )





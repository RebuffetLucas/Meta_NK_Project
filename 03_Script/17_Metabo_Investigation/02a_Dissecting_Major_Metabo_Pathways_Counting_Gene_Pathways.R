
PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds") #Object with all the genes


#Attention, à refaire avec tous les gènes inclus!
library(devtools)

#install_github("guokai8/scGSVA")
set.seed(123)   
library(scGSVA)   

#Load res object
res= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/res_BasedOn23Kgenes.rds") #based on the object with all the genes

#Build Annotation
hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")

Heatmap(res,group_by="FirstClust", average= TRUE)
dotPlot(res,features="Wnt.signaling.pathway",group_by="FirstClust")


str(hsko)


#Investigate Metabo

METABO_PATHWAY_INTEREST = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Pahway_Metabo_To_Investigate_Metabo_Interest2.txt", header= FALSE, quote= "\"" )
METABO_PATHWAY_INTEREST_Weird_Name = METABO_PATHWAY_INTEREST[,"V1"]


#Convert weird names into Normal names
Annotations_Conversion_DF = data.frame(unique(res@annot[["Annot"]]))
Annotations_Conversion_DF= na.omit(Annotations_Conversion_DF)
Annotations_Conversion_DF$Weird_Name=colnames(res@gsva)

colnames(Annotations_Conversion_DF) = c("Normal_Name", "Weird_Name")

Annotations_Conversion_DF %>%  filter(Weird_Name %in% METABO_PATHWAY_INTEREST_Weird_Name) -> Annotations_Of_Interest

res@annot %>%  filter(Annot %in% Annotations_Of_Interest$Normal_Name) -> Genes_Annotated

Genes_Dataset = rownames(PBMC)

# Create a dataframe summary of the pathways of interest
Summary_Pathways <- data.frame(
  Annot = unique(Genes_Annotated$Annot),
  GeneID_List = sapply(unique(Genes_Annotated$Annot), function(x) toString(Genes_Annotated$GeneID[Genes_Annotated$Annot == x])),
  Number_Genes_Pathway = sapply(unique(Genes_Annotated$Annot), function(x) length(Genes_Annotated$GeneID[Genes_Annotated$Annot == x]))
)
 

Summary_Pathways

write.csv(Summary_Pathways, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Summary_Pathways.csv", row.names=FALSE)

#Identify which enzymes are differentially expressed from one cluster to the other

PBMC= SetIdent(PBMC, value= "FirstClust")

#Extract all the metabo genes
  # Extract unique GeneID lists
unique_geneid_lists <- unique(Summary_Pathways$GeneID_List)
  # Merge all unique GeneID lists into one single list
Genes_To_Test <- unlist(strsplit(unique_geneid_lists, ", "))
Genes_To_Test2 <- strsplit(unique_geneid_lists, ", ") # Can be used to look in heatmap etc...


#test which one are differentially expressed DOES NOT WORK, 
#All_Markers = FindAllMarkers(PBMC , features= Genes_To_Test, only.pos = FALSE, verbose = TRUE, min.pct =  0 , logfc.threshold = 0.5 )
#All_Markers = FindAllMarkers(PBMC , features= Genes_To_Test, only.pos = FALSE, method= FINDMARKERS_METHOD , min.pct =  0.01 , logfc.threshold = 0.1 , verbose = TRUE)

#Have a look at the level of expression

VlnPlot(PBMC, features= Genes_To_Test2[[7]])

# Make a Complex Heatmap divided by Pathway only with genes present in more than 5% cells








#Subset values of Interest
df_Metabo = res@gsva[, METABO_PATHWAY_INTEREST$V1]



Heatmap(res , features= METABO_PATHWAY_INTEREST$V1, group_by = "FirstClust", average = TRUE)

palette_pheatmap =  c(   NK3 ="#ABA300", NK1C =  "#F8766D"   , NK1A = "#0CB702"  , NK1B  =  "#00BFC4"   ,   NK2=    "#ED68ED"   , NKint= "#8494FF"  )

palette_pheatmap2 =  c(    "#F8766D"   ,"#0CB702"  ,  "#00BFC4"     ,  "#8494FF" ,   "#ED68ED" , "#ABA300" )

ORDER_LEGEND= c("NK1C", "NKint", "NK1B", "NK2", "NK1A", "NK3" )
ORDER_LEGEND= c( NK1C, NKint, NK1B, NK2, NK1A, NK3 )

p = Heatmap(res , features= METABO_PATHWAY_INTEREST$V1, group_by = "FirstClust", average = TRUE, scale = "row", fontsize_row = 10, fontsize_col = 12 , annotation_colors =  list(palette_pheatmap2), angle_col = 45) 


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Pathway_Metabo_Heatmap.png", width = 20, height = 15,  units = "cm", res=600 )
p
dev.off()





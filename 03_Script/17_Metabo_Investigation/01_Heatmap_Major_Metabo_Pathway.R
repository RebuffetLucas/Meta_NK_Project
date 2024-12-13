
PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds") #Object with all the genes


#Attention, à refaire avec tous les gènes inclus!
library(devtools)

install_github("guokai8/scGSVA")
set.seed(123)   
library(scGSVA)   

####### Run only the first time #########
hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
res<-scgsva(PBMC,hsko,method="ssgsea") ## or use UCell
Heatmap(res,group_by="FirstClust")
dotPlot(res,features="Wnt.signaling.pathway",group_by="groups")
##########





############## All 33 Metabo of INTEREST ###################
METABO_PATHWAY_INTEREST = read.csv("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Pahway_Metabo_To_Investigate_Metabo_Interest2.txt", header= FALSE, quote= "\"" ) #30 pathways = start

############## All 33 Metabo of INTEREST ###################
METABO_PATHWAY_INTEREST = read.csv("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Pahway_Metabo_To_Investigate_Pathway_Reliable_Only.txt", header= FALSE, quote= "\"" ) #19 pathways = reliable



#saveRDS(res, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/res_BasedOn23Kgenes.rds")

#res_pas_ouf = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/res_First_Try.rds") #Based on the obkect with 13K genes
res= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/res_BasedOn23Kgenes.rds") #based on the object with all the genes


#write.table (colnames(res@gsva), file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/pathways.txt")


res@gsva[1:10, 1:10]
res_pas_ouf@gsva[1:10, 1:10]



str(res@gsva)
head(res@gsva)

colnames(res@gsva)
grepl("metabo" , names(res@gsva))

res@gsva


help("Heatmap")

METABO_PATHWAY_INTEREST = METABO_PATHWAY_INTEREST[,"V1"]

#Subset values of Interest
df_Metabo = res@gsva[, METABO_PATHWAY_INTEREST]

Heatmap(res , features= METABO_PATHWAY_INTEREST, group_by = "FirstClust", average = TRUE)

p = Heatmap(res , features= METABO_PATHWAY_INTEREST, group_by = "FirstClust", average = TRUE)




palette_pheatmap =  c(   NK3 ="#ABA300", NK1C =  "#F8766D"   , NK1A = "#0CB702"  , NK1B  =  "#00BFC4"   ,   NK2=    "#ED68ED"   , NKint= "#8494FF"  )

palette_pheatmap2 =  c(    "#F8766D"   ,"#0CB702"  ,  "#00BFC4"     ,  "#8494FF" ,   "#ED68ED" , "#ABA300" )

ORDER_LEGEND= c("NK1C", "NKint", "NK1B", "NK2", "NK1A", "NK3" )

p = Heatmap(res , features= METABO_PATHWAY_INTEREST, group_by = "FirstClust", average = TRUE, scale = "row", fontsize_row = 10, fontsize_col = 12 , angle_col = 45) 



png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Pathway_Metabo_Heatmap_33Pathways_Interest.png", width = 20, height = 15,  units = "cm", res=600 )
p
dev.off()


png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Pathway_Metabo_Heatmap_19Reliable_PathwaysOfInterest.png", width = 20, height = 15,  units = "cm", res=600 )
p
dev.off()

#Fig3d
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig3d.pdf",
  plot = p,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 20,
  height = 14,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



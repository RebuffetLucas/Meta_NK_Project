#Apply scores of Dim, Bright and Adapt to the initial UMAP

data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#All

PBMC= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

DimPlot(PBMC, label = TRUE, cols= palette)

PBMC2= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes.rds")
PBMC2@meta.data = PBMC@meta.data

PBMC2@reductions[["UMAP"]] = PBMC@reductions[["UMAP"]]

PBMC_Norm= as.SingleCellExperiment(PBMC2)
PBMC_Norm = multiBatchNorm(PBMC_Norm,  batch = PBMC_Norm$orig.ident)
PBMC2= as.Seurat(PBMC_Norm)

PBMC2$FirstClust = factor(PBMC$FirstClust, levels= ORDER_CLUST_LEGEND)

DimPlot(PBMC2, group.by = "FirstClust", cols = palette)


#PBMC New Names + All Genes
PBMC2= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

PBMC2 = SetIdent(PBMC2, value="FirstClust")
FILTER_pvalue_THR = 0.05

rownames(PBMC2)[grepl("KIR" , rownames(PBMC2))]


Human_KIR= c("KIR2DL3" ,"KIR2DL1", "KIR3DL1" ,"KIR3DL2", "KIR2DL4", "KIR3DL3")


VlnPlot(PBMC2, features= Human_KIR, group.by = "FirstClust")





#Do the ILCP scoring from Colonna

#Try ILCP from Colonna:

Markers_Fig1a = read.csv("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/markers_Fig1a.csv")


Markers_Fig1a %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_Fig1a


top_Fig1a %>% filter(cluster== "ILCP") -> Top_ILCP

PBMC2 = AddModuleScore( PBMC2 , features = list(Top_ILCP$gene), pool= NULL ,name= "ILCP_Colonna" , seed=19) 



#Import signatures


FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Bhandoola_V2_raw.xlsx")

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



#Monitored markers
MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)
MONITORED_Markers  = MONITORED_Markers$fm2_ENKP_ILCP

MONITORED_Markers %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::arrange(p_val_adj) %>%
  filter(p_val_adj<FILTER_pvalue_THR)  -> top40_ENKP



top40_ENKP$gene


#Convert into human genes
#Sabrina's function to convert to human genes

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  mouse_attributes <- c("external_gene_name","external_gene_source","hsapiens_homolog_associated_gene_name","hsapiens_homolog_perc_id",
                        "hsapiens_homolog_goc_score", "hsapiens_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "hgnc_symbol",
                  filters = "hgnc_symbol", values = x, mart = ensembl.human,
                  attributesL = mouse_attributes, martL = ensembl.mouse,uniqueRows = TRUE)
  mousex <- genesV2 %>%
    filter(Gene.name != "") %>%
    group_by(HGNC.symbol) %>%
    slice_max(Human.Gene.order.conservation.score) %>%
    filter(Human.orthology.confidence..0.low..1.high. == 1) %>%
    dplyr::select(HGNC.symbol,Gene.name) %>%
    dplyr::rename(human.gene.name = HGNC.symbol, mouse.gene.name = Gene.name)
  
  return(mousex)
}


convertMouseGeneList <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  human_attributes <- c("external_gene_name","external_gene_source","mmusculus_homolog_associated_gene_name","mmusculus_homolog_perc_id",
                        "mmusculus_homolog_goc_score", "mmusculus_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "mgi_symbol",
                  filters = "mgi_symbol", values = x, mart = ensembl.mouse,
                  attributesL = human_attributes, martL = ensembl.human,uniqueRows = TRUE)
  
  humanx <- genesV2 %>%
    filter(Gene.name != "") %>%
    group_by(MGI.symbol) %>%
    slice_max(Mouse.Gene.order.conservation.score) %>%
    filter(Mouse.orthology.confidence..0.low..1.high. == 1) %>%
    dplyr::select(MGI.symbol,Gene.name) %>%
    dplyr::rename(mouse.gene.name = MGI.symbol, human.gene.name = Gene.name)
  
  return(humanx)
}

  #Convert into human
musGenes_ENKP <- top40_ENKP$gene
Human_ENKP = convertMouseGeneList(musGenes_ENKP)


top40_ENKP <- top40_ENKP %>%
  left_join(Human_ENKP, by = c("gene" = "mouse.gene.name"))

top40_ENKP[top40_ENKP$gene=="Ccl5","human.gene.name"]= "CCL5"

top40_ENKP[top40_ENKP$gene=="Klra8","human.gene.name"]= Human_KIR[1]
top40_ENKP[top40_ENKP$gene=="Klri2","human.gene.name"]= Human_KIR[2]
top40_ENKP[top40_ENKP$gene=="Klra9","human.gene.name"]= Human_KIR[3]
top40_ENKP[top40_ENKP$gene=="Klra3","human.gene.name"]= Human_KIR[4]
top40_ENKP[top40_ENKP$gene=="Klra7","human.gene.name"]= Human_KIR[5]
top40_ENKP[top40_ENKP$gene=="Klra4","human.gene.name"]= Human_KIR[6]


top40_ENKP %>%
  filter(!is.na(human.gene.name))%>%
  slice_head(n = 20)  -> top20_ENKP_human

#write.csv(top20_ENKP_human , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/18_Table_Sup/Top20_ENKP.csv", row.names=FALSE) #Save for Data sup


#write.csv(top20_ILCP_human , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/18_Table_Sup/Top20_ILCP.csv", row.names=FALSE) #Save for Data sup


top20_ENKP_human2 = list(intersect(rownames(PBMC2), top20_ENKP_human$human.gene.name) )


length(unlist(top20_ENKP_human2))


DotPlot(PBMC2 , features= top20_ENKP_human$human.gene.name, group.by = "FirstClust")

DotPlot(PBMC2 , features= unlist(top20_ENKP_human2), group.by = "FirstClust")

#Save the signatures used for scoring
data_export = data.frame(top20_ILCP_human$gene , top20_ILCP_human$human.gene.name ,top20_ENKP_human$gene  ,top20_ENKP_human$human.gene.name)
names(data_export)= c("top20_ILCP_mouse","top20_ILCP_human", "top20_ENKP_human", "top20_ENKP_mouse" )

write.csv(data_export, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/Progenitors_SCoring/ENKP_ILCP_genes.csv", row.names=FALSE)



#Keep only the one included in rownames
#MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(PBMC)))
#MONITORED_Markers= lapply(MONITORED_Markers, FUN=function(x) head(x, n=40))

#Add the module scores

PBMC2 = AddModuleScore(PBMC2, features = top20_ENKP_human2 , name = "ENKP")

VlnPlot(PBMC2, features = "ILCP_Colonna1")
VlnPlot(PBMC2, features = "ENKP1")



#FeaturePlot ENKP and ILCP



#ENKP
p5 = FeaturePlot(PBMC2, features = "ENKP1" , pt.size= 0.6, max.cutoff = 1)   &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) ) & NoAxes()

FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/pval_", FILTER_pvalue_THR,"_ENKP","FeaturePlot","_Reds",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p5
dev.off()



p6 = FeaturePlot(PBMC2, features = "ENKP1" , pt.size= 0.6, max.cutoff = 1)   &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/pval_", FILTER_pvalue_THR,"_ENKP","FeaturePlot","RdBu",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p6
dev.off()



#Fig5a
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig5a.pdf",
  plot = p6,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)




p15 = FeaturePlot(PBMC2, features = "ILCP_Colonna1" , pt.size= 0.6, max.cutoff = 1)   &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) ) & NoAxes()

FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/pval_", FILTER_pvalue_THR,"_ILCP_Colonna1","FeaturePlot","_Reds",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p5
dev.off()



p16 = FeaturePlot(PBMC2, features = "ILCP_Colonna1" , pt.size= 0.6, max.cutoff = 1)   &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/pval_", FILTER_pvalue_THR,"_ILCP_Colonna1","FeaturePlot","RdBu",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p16
dev.off()


#Fig5c
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig5c.pdf",
  plot = p16,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



# ViolInPlot ENKP and ILCP
p1 = VlnPlot(PBMC2, features = "ILCP_Colonna1", group.by = "FirstClust", pt.size = 0, cols = palette, y.max = 1.1) & stat_summary(fun.data=data_summary,color="black")

FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/pval_", FILTER_pvalue_THR,"_ILCP_Colonna","VlnPlot",".png")

png(file=FILE, width = 20, height = 15,  units = "cm", res=600 )
p1
dev.off()



#Fig5d
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig5d.pdf",
  plot = p1,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)







p2 = VlnPlot(PBMC2, features = "ENKP1", group.by = "FirstClust", pt.size=0, cols = palette) & stat_summary(fun.data=data_summary,color="black")


FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/pval_", FILTER_pvalue_THR,"_ENKP","VlnPlot",".png")

png(file=FILE, width = 20, height = 15,  units = "cm", res=600 )
p2
dev.off()


#Fig5b
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig5b.pdf",
  plot = p2,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



p11 = RidgePlot(PBMC2, features = "ENKP1" , group.by= "FirstClust", cols= palette, sort = TRUE)  

FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/","ENKP","RidgePlot",".png")

png(file=FILE, width = 30, height = 20,  units = "cm", res=600 )
p11
dev.off()


p12 = RidgePlot(PBMC2, features = "ILCP_Colonna1" , group.by= "FirstClust", cols= palette, sort = TRUE)  

FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/","ILCP_Colonna1","RidgePlot",".png")

png(file=FILE, width = 30, height = 20,  units = "cm", res=600 )
p12
dev.off()








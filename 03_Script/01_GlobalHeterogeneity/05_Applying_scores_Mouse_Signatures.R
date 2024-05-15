#Apply scores of Dim, Bright and Adapt to the initial UMAP

#All
#read Object
PBMC2= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

PBMC2 = SetIdent(PBMC2, value="FirstClust")
FILTER_pvalue_THR = 0.05

rownames(PBMC2)[grepl("KIR" , rownames(PBMC2))]


Human_KIR= c("KIR2DL3" ,"KIR2DL1", "KIR3DL1" ,"KIR3DL2", "KIR2DL4", "KIR3DL3")


#Import signatures


FILE_SIGNATURES_PATH = file.path( MOUSE_SIGNATURE_xlsx_path)

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

MONITORED_Markers %>%
  dplyr::filter(avg_log2FC < 0) %>%
  dplyr::arrange(p_val_adj) %>%
  filter(p_val_adj<FILTER_pvalue_THR) -> top40_ILCP


top40_ENKP$gene
top40_ILCP$gene

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

#Manually add some genes unasigned
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



musGenes_ILCP <- top40_ILCP$gene
Human_ILCP = convertMouseGeneList(musGenes_ILCP)

VlnPlot(PBMC2, feature= "CD7", group.by = "FirstClust")

top40_ILCP <- top40_ILCP %>%
  left_join(Human_ILCP, by = c("gene" = "mouse.gene.name"))

top40_ILCP[top40_ILCP$gene=="Cd7","human.gene.name"]= "CD7"
top40_ILCP[top40_ILCP$gene=="Klrb1f","human.gene.name"]= "KLRB1"

top40_ILCP %>%
  filter(!is.na(human.gene.name))%>%
  slice_head(n = 20)  -> top20_ILCP_human



top20_ILCP_human2  =  list(  intersect(rownames(PBMC2), top20_ILCP_human$human.gene.name) )
top20_ENKP_human2 = list(intersect(rownames(PBMC2), top20_ENKP_human$human.gene.name) )

length(unlist(top20_ILCP_human2))
length(unlist(top20_ENKP_human2))

DotPlot(PBMC2 , features= top20_ILCP_human$human.gene.name, group.by = "FirstClust")
DotPlot(PBMC2 , features= top20_ENKP_human$human.gene.name, group.by = "FirstClust")

DotPlot(PBMC2 , features= unlist(top20_ILCP_human2), group.by = "FirstClust")
DotPlot(PBMC2 , features= unlist(top20_ENKP_human2), group.by = "FirstClust")

#Save the signatures used for scoring
data_export = DataFrame(top20_ILCP_human$gene , top20_ILCP_human$human.gene.name ,top20_ENKP_human$gene  ,top20_ENKP_human$human.gene.name)
names(data_export)= c("top20_ILCP_mouse","top20_ILCP_human", "top20_ENKP_human", "top20_ENKP_mouse" )

write.csv(data_export, PATH_csv_ENKP_ILCP_Genes, row.names=FALSE)



#Keep only the one included in rownames
#MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(PBMC)))
#MONITORED_Markers= lapply(MONITORED_Markers, FUN=function(x) head(x, n=40))

#Add the module scores

PBMC2 = AddModuleScore(PBMC2, features = top20_ILCP_human2 , name = "ILCP") 
PBMC2 = AddModuleScore(PBMC2, features = top20_ENKP_human2 , name = "ENKP")

#Dataviz
VlnPlot(PBMC2, features = "ILCP1")
VlnPlot(PBMC2, features = "ENKP1")


#Save FeaturePlot ENKP and ILCP

p3 = FeaturePlot(PBMC2, features = "ILCP1" , pt.size= 0.6)   &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) ) & NoAxes()

FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ILCP","FeaturePlot","_Reds",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p3
dev.off()


p4 = FeaturePlot(PBMC2, features = "ILCP1" , pt.size= 0.6)   &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )), midpoint = 0 ) & NoAxes()


FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ILCP","FeaturePlot","_RdBu",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p4
dev.off()


p12 = RidgePlot(PBMC2, features = "ILCP1" , group.by= "FirstClust", cols = palette)  

FILE= paste0("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/","ILCP","RidgePlot",".png")

png(file=FILE, width = 30, height = 20,  units = "cm", res=600 )
p12
dev.off()




#ENKP
p5 = FeaturePlot(PBMC2, features = "ENKP1" , pt.size= 0.6)   &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) ) & NoAxes()

FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ENKP","FeaturePlot","_Reds",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p5
dev.off()

p6 = FeaturePlot(PBMC2, features = "ENKP1" , pt.size= 0.6)   &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()

FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ENKP","FeaturePlot","_RdBu",".png")

png(file=FILE, width = 22, height = 18,  units = "cm", res=600 )
p6
dev.off()




# ViolinPlot ENKP and ILCP
p1 = VlnPlot(PBMC2, features = "ILCP1", group.by = "FirstClust", pt.size = 0, cols = palette)

FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ILCP","VlnPlot",".png")

png(file=FILE, width = 20, height = 20,  units = "cm", res=600 )
p1
dev.off()


p2 = VlnPlot(PBMC2, features = "ENKP1", group.by = "FirstClust", pt.size=0, cols = palette)


FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ENKP","VlnPlot",".png")

png(file=FILE, width = 20, height = 20,  units = "cm", res=600 )
p2
dev.off()

p11 = RidgePlot(PBMC2, features = "ENKP1" , group.by= "FirstClust", cols= palette)  

FILE= paste0(PATH_FIGURES, FILTER_pvalue_THR,"_ENKP","RidgePlot",".png")

png(file=FILE, width = 30, height = 20,  units = "cm", res=600 )
p11
dev.off()



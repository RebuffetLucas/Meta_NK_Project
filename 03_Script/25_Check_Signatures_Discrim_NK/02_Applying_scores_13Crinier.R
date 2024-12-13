#Import signatures
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Crinier_13_NKgenes.xlsx") #13 gèes Crinier 2018

#def function

# +- std
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}



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




MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)

MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(PBMC)))
MONITORED_Markers= lapply(MONITORED_Markers, FUN=function(x) head(x, n=40))

#Add the module scores
for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

#Convert to factors:
PBMC$celltype.l1 = as.factor(PBMC$celltype.l1)
PBMC$celltype.l2 = as.factor(PBMC$celltype.l2)
PBMC$celltype.l3 = as.factor(PBMC$celltype.l3)

# At level 1
PBMC = SetIdent(PBMC,  value= "celltype.l1")
p1 = VlnPlot(PBMC,features = "Crinier_13_hNKGenes_20181", pt.size = 0, sort = "decreasing" ) & ggtitle("Score 13 NK Genes")

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_13NKgenes_OnCD45pos_l1.png"), width = 25, height = 15,  units = "cm", res=600 )
p1
dev.off()

  #Increasing
p1bis = VlnPlot(PBMC,features = "Crinier_13_hNKGenes_20181", pt.size = 0, sort = "decreasing") & ggtitle("Score 13 NK Genes") & stat_summary(fun.data=data_summary,color="black")

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_13NKgenes_OnCD45pos_l1_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
p1bis
dev.off()



# At level 2
PBMC = SetIdent(PBMC,  value= "celltype.l2")

p2 = VlnPlot(PBMC,features = "Crinier_13_hNKGenes_20181", pt.size = 0, sort = "decreasing") & ggtitle("Score 13 NK Genes")

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_13NKgenes_OnCD45pos_l2.png"), width = 30, height = 15,  units = "cm", res=600 )
p2
dev.off()


p2bis = VlnPlot(PBMC,features = "Crinier_13_hNKGenes_20181", pt.size = 0, sort= "decreasing") & ggtitle("Score 13 NK Genes") & stat_summary(fun.data=data_summary,color="black")

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_13NKgenes_OnCD45pos_l2_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
p2bis
dev.off()




# At level  3
PBMC = SetIdent(PBMC,  value= "celltype.l3")

p3 = VlnPlot(PBMC,features = "Crinier_13_hNKGenes_20181", pt.size = 0,  sort = "decreasing") & ggtitle("Score 13 NK Genes")

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_13NKgenes_OnCD45pos_l3.png"), width = 60, height = 15,  units = "cm", res=600 )
p3
dev.off()


p3bis = VlnPlot(PBMC,features = "Crinier_13_hNKGenes_20181", pt.size = 0, sort= "decreasing") & ggtitle("Score 13 NK Genes") & stat_summary(fun.data=data_summary,color="black")

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_13NKgenes_OnCD45pos_l3_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
p3bis
dev.off()




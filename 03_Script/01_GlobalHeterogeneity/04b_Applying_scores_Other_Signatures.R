#Apply scores of Dim, Bright and Adapt to the initial UMAP

#All

data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#Ra=ead Object and set ident
PBMC= readRDS(PATH_CURATED_Object)
PBMC= SetIdent(PBMC, value= "FirstClust")

#dataviz
DimPlot(PBMC, label = TRUE, cols = palette)



#Import signatures
FILE_SIGNATURES_PATH = PATH_xlsx_Signature_of_interest #Path to the xlsx of signatures


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



#Turn it into a list suitable for AddModuscore
MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)

MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(PBMC)))
MONITORED_Markers= lapply(MONITORED_Markers, FUN=function(x) head(x, n=NUMBER_MARKERS_SCORING))

#Add the module scores
for (i in names(MONITORED_Markers)){
  PBMC = AddModuleScore(PBMC, features = list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}





##### Vizualization #####


#FeaturePlot

p1 = FeaturePlot(PBMC, features= SIGNATURE_NAME, pt.size = 0.6) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
png(file=paste0(PATH_OUTPUT , SIGNATURE_NAME, "FeaturePLOT_MODULE_score.png"), width = 18, height = 18,  units = "cm", res=600 )
p1
dev.off()




#Violin Plot
p2 = VlnPlot(PBMC , features = paste0(names(MONITORED_Markers),1), group.by= "FirstClust" , pt.size = 0, cols = palette)  & stat_summary(fun.data=data_summary,color="black")
png(file=paste0(PATH_OUTPUT , SIGNATURE_NAME, "Vln_Plot_MODULE_score.png"), width = 18, height = 18,  units = "cm", res=600 )
p2
dev.off()



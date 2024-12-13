#Running the Label Transfer process on Tang et.al.

knitr::opts_chunk$set(fig.align="center",eval.after="fig.cap",
                      fig.width=12,fig.height=5,
                      warning=FALSE, 
                      results="asis",
                      message=FALSE)

###############################################
############## Def a few functions ############
###############################################


## scale rows ####
# as scale="row" does not seem to work well in pheatmap... 
#  https://github.com/raivokolde/pheatmap/blob/44a44f7e640c73867f40261691c16dfa4da34fa8/R/pheatmap.r
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

## Object name --> caracter string
name2string <- function(v1) {
  deparse(substitute(v1))
}


# Vector to sentence
VecToEnumeration<-function(vec) {
  
  if (length(vec)>2) {
    vec1<-vec[1:(length(vec)-1)]
    en1<-paste0(vec1,collapse = ", ")
    en2<-paste0(en1, " and ", vec[length(vec)])
  }
  if (length(vec)==2) {
    en2<-paste0(vec[1], " and ", vec[2])
  }
  if (length(vec)==1) {
    en2<-vec
  }
  return(en2)
}

# darken colors https://gist.github.com/Jfortin1/72ef064469d1703c6b30
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

FeaturePlot2<-function(object,feature,reduc="umap",dims=1:2) {
  exp <- GetAssayData(object, "data")[feature,]
  coor<-(object@reductions[[reduc]])@cell.embeddings[,dims]
  coor<-coor[order(exp),]
  exp<-sort(exp)  
  
  if (length(exp)<10^5) {
    plot(coor,col=exp2col(exp),pch=19,cex=0.2, main=feature)
  } else {
    data<-data.frame(exp=exp,coor)
    colNames<-colnames(data)
    p<-ggplot(data = data)+geom_scattermore( 
      mapping = aes_string(x = colNames[2], 
                           y = colNames[3], 
                           color = exp),
      pointsize = 0.8, pixels = c(1000, 1000)) + 
      ggtitle(g)+ 
      scale_color_gradientn(colors=c("royalblue","springgreen","yellow","red")) + 
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 
    
    print(p)
    
  }
  
}

exp2col<-function(value,breaks=1000,
                  cols=c("royalblue","springgreen","yellow","red")){
  colPal <- colorRampPalette(cols)
  marker.color <- colPal(breaks)[as.numeric(cut(value, breaks = breaks))]
  return(marker.color)
}

# dotPlot2 ####
DotPlot2=function(object,group.by = NULL, features,scale=TRUE,
                  dot.min=0,dot.scale = 6) {
  p=DotPlot(object=object,group.by=group.by,
            features=features, scale=scale,
            dot.min=dot.min,dot.scale = dot.scale )+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("") + ylab("Cluster")+
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust=1,
                                     size=8))+
    scale_colour_gradientn(colours =  c("#1c0142", "#501dcf", "#5f9cc9",
                                        "#72dbae", "#92c77d", "#FFDB00",
                                        "#FF0000"))
  return(p)
}


## extract markers #
ExtractMarkers=function(DEresults) {
  markers=c()
  for (cl in sort(unique(DEresults$cluster))) {
    print(cl)
    markers_overCl_cl=DEresults[
      DEresults$avg_log2FC>0&
        DEresults$p_val_adj<=0.05&
        DEresults$cluster==cl,
    ]
    
    markers[[cl]]=markers_overCl_cl[,"gene"]
    
  }
  return(markers)
}


# removing duplicated markers
RemoveDup=function(features_list) {
  features=unlist(features_list)
  dup=features[duplicated(features)]
  
  for (cl in names(features_list)) {
    genes=features_list[[cl]]
    genes=genes[!genes%in%dup]
    features_list[[cl]]=genes
  }
  return(features_list)
  
}


# extract top markers #
extractTopMarkers=function(markers,n=10) {
  markers %>%
    filter(avg_log2FC>0) %>%
    filter(p_val_adj<=0.05) %>%
    group_by(cluster) %>%
    data.frame()  -> markers_pos_markers
  markers_pos_markers$cluster=as.character(markers_pos_markers$cluster)
  
  top_markers=c()
  for (cl in unique(markers_pos_markers$cluster)) {
    markers_cl=markers_pos_markers[
      markers_pos_markers$cluster==cl,
    ]
    markers_cl=markers_cl[order(markers_cl$p_val),]
    markers_cl_gene=markers_cl$gene[1:n]
    markers_cl_gene=markers_cl_gene[!is.na(markers_cl_gene)]
    top_markers[[cl]]=  markers_cl$gene[1:n]
  }
  #top_markers=unique(top_markers)
  
  return(top_markers)
}


#######################################
############ Import Libraries #########
#######################################

library(DESeq2)
library(edgeR)
library(ade4)
library(factoextra)
library(kableExtra)
library(RColorBrewer)
library(vsn)
library(ggVennDiagram)
library(wesanderson)
library(DescTools)
library(pheatmap)
library(WGCNA)
library(readxl)
library(gage)
library(sva)
library(gageData)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(ggrepel)
library(openxlsx) 
library(ggrastr)
library(viridis)
library(kableExtra)
library(msigdbr)
library(pathview)
data(gene.idtype.list)
library(biomaRt)
library(SBGNview)
library(limma)
library(filesstrings)
#library(GSVA)
#library(fgsea)
library(DT)
#library(ComplexUpset)
library(grid)
library("ggVennDiagram")
#library(WebGestaltR)
library(biomaRt)
require(VennDiagram)
library("extrafont")
library("ggvenn")
library(Seurat)
library(scattermore)
library(UniprotR)
library(reticulate)
library(harmony)


##########################################
############# Harmony Method #############
##########################################

#Load data
PBMC_Meta_NK = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")
DimPlot(PBMC_Meta_NK)
PBMC_Meta_NK

genes_to_keep = intersect(rownames(PBMC_Meta_NK), rownames(PBMC_Blood_All))

PBMC_All = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5seurat")

PBMC_Blood_All = subset(PBMC_All, subset = meta_tissue ==  "Blood" & meta_platform == "10X" )



#Merge data
data_all =  merge(PBMC_Blood_All, y = PBMC_Meta_NK, add.cell.ids = c("Tang", "Meta"))
data_all = subset(data_all, features= genes_to_keep)
data_all

#Meta Data
# metadata metaNK ####
metadata_NK=PBMC_Meta_NK@meta.data[,c("Dataset","orig.ident","FirstClust","SecondClust")]
metadata_NK$Project="MetaNK"
metadata_NK$meta_histology = as.factor("Healthy donor")
rownames(metadata_NK) = paste0("Meta_",rownames(metadata_NK))
# metadata Tang ####
metadata_Tang=PBMC_Blood_All@meta.data[,c("meta_patientID","meta_histology","celltype","Majortype")]
colnames(metadata_Tang)[1]="orig.ident"
metadata_Tang$FirstClust=metadata_Tang$Majortype
metadata_Tang$SecondClust=metadata_Tang$celltype
metadata_Tang$Dataset=metadata_Tang$orig.ident
metadata_Tang$Project="Tang"
rownames(metadata_Tang) = paste0("Tang_",rownames(metadata_Tang))
metadata_Tang=metadata_Tang[,
                          c("Dataset","orig.ident",
                            "FirstClust","SecondClust","Project", "meta_histology")]



# Combining metadata ####
metadata_all=rbind(metadata_NK,metadata_Tang)
metadata_all=metadata_all[colnames(data_all),]


data_all@meta.data = metadata_all


data_all  = NormalizeData(data_all)


#Variable Features:
VF1 = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VariableFeatures.rds")

intersect(VF1, genes_to_keep)

intersect(VariableFeatures(data_all) , intersect(VF1, genes_to_keep) )
VariableFeatures(data_all)=intersect(VF1, genes_to_keep)
#data_all = FindVariableFeatures(data_all, selection.method = "vst", nfeatures = 2000)

data_all = ScaleData(data_all, do.scale = DATA_SCALE, do.center = DATA_CENTER)

data_all= RunPCA(data_all, npcs = 20 )

DimPlot(data_all, reduction= "pca", group.by = "Project")

data_all=data_all %>%                  
  RunHarmony(c("Dataset"), plot_convergence = TRUE, reduction = "pca" ,  verbose = TRUE , dims.use= 1:20 ) 

data_all <- RunUMAP(data_all, reduction = "harmony", dims= 1:20, reduction.name = "UMAP", reduction.key = "UMAP")

## umap in the good sign
umap_coord=data_all@reductions$UMAP@cell.embeddings
median_x0=aggregate(umap_coord[,1],list(data_all$FirstClust),FUN=median)
median_x=median_x0[,2]%>% setNames(median_x0[,1])
if (median_x["NK2A"]>median_x["NK1C"]) {
  umap_coord[,1]=(-umap_coord[,1])
}
median_y0=aggregate(umap_coord[,2],list(data_all$FirstClust),FUN=median)
median_y=median_y0[,2]%>% setNames(median_y0[,1])
if (median_y["NK2A"]<median_y["NK1B"]) {
  umap_coord[,2]=(-umap_coord[,2])
}
data_all[["umap"]] <- CreateDimReducObject(embeddings = umap_coord, key = "UMAP_", assay = DefaultAssay(data_all))



DimPlot(data_all, label= TRUE, group.by = "FirstClust", split.by = "Project", reduction = "umap")
DimPlot(data_all, label= TRUE, group.by = "SecondClust", split.by = "Project", reduction = "umap")
DimPlot(data_all, label= TRUE, group.by = "SecondClust", split.by = "Project")



#data_all <- RunUMAP(data_all,  reduction.name = "UMAP", reduction.key = "UMAP")

data_all <- FindNeighbors(data_all, reduction = "harmony", dims= 1:20)
data_all <- FindClusters(data_all, resolution = 0.8, random.seed = SEED )

#saveRDS(data_all, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_With_TangPaper/Integrated_Seurat_Object/IntegratedMetaV2_Tang.rds")
DimPlot(data_all, label= TRUE, split.by =  "Project")



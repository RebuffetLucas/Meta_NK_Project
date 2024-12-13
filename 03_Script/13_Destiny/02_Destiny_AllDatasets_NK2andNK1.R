
#Install Libraries
#Works
#install.packages("sinaplot")
#BiocManager::install("destiny")

#install.packages("rgl")
#install.packages("magick")
#BiocManager::install("genefilter")

#Does not work
#remotes::install_github("Japrin/sscVis")
#devtools::install_github("Japrin/sscVis")

ROOT_To_CHOOSE = 1


#Library Import

library(Seurat)
library(dplyr)
#library(sscVis)
library(harmony)
#library(UCell)
library(patchwork)
#library(Matrix.utils)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(sinaplot)
library(destiny)
library(rgl)
library(magick)
library(genefilter)
library(SeuratWrappers)
library("scatterplot3d") 
library(batchelor)


#Loading and preparing the data ###

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

object= PBMC
DimPlot(PBMC)
var.features = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VariableFeatures.rds")

set.seed(SEED)

col_firstCl=c(NK1A = "#0CB702", 
              NK1B = "#00BFC4", NK1C = "#F8766D",
              NK2A = "#ED68ED", NK2B = "#8494FF", 
              NK3 = "#ABA300")



log.norm.data=object@assays$RNA@data
log.norm.data.var<-log.norm.data[var.features,]
data<-t(as.matrix(log.norm.data.var))


metadata=object@meta.data
metadata$Dataset2=as.character(metadata$Dataset)
metadata$Dataset2[metadata$orig.ident%in%c("CMVneg1_donorA",
                                           "CMVneg2_donorB")]="Dataset4_neg"
metadata$Dataset2[metadata$orig.ident%in%c("CMVpos2_donorC",
                                           "CMVpos3_donorD", 
                                           "CMVpos4_donorE")]="Dataset4_pos"

#Diffusion map (FastMNN)

#Subset
#PBMC2  =subset(PBMC, subset= Dataset == "Dataset4")
PBMC2 = subset(PBMC, subset = FirstClust== "NK3", invert= TRUE)

DimPlot(PBMC2)

VariableFeatures(PBMC2) = var.features

object_d = PBMC2

object_d2 <- RunFastMNN(object.list = SplitObject(object_d, 
                                                  split.by = "orig.ident"),
                        features=var.features,
                        assay = "RNA")

matrix=t(object_d2@assays$mnn.reconstructed@data) %>% as.matrix()


#dm= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/DiffusionMap_FastMNN_AllDatasets_NK1and2.rds")


dm<-DiffusionMap(data=matrix,
                 censor_val = 30,
                 censor_range = c(30, 40),
                 verbose=TRUE
)

plot(dm)

#Save

#saveRDS(dm , "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/DiffusionMap_FastMNN_AllDatasets_NK1and2.rds")

dm_coor=dm@eigenvectors


#Root with automatic choice:
dpt_d <- DPT(dm)





plot(dpt_d, root = 1)
plot(dpt_d, root = 2)
plot(dpt_d, root = 3)


root_list <- which(apply(dpt_d@tips, 1, function(x) all(x == TRUE, na.rm=TRUE)))


root  = root_list[ROOT_To_CHOOSE]


#Root with NK2A as initial
group=factor(metadata[rownames(dm_coor),"FirstClust"], levels=c("NK2A","NK2B","NK1A","NK1B","NK1C"))

y=aggregate(dm_coor[,1],by=list(group),median)
y2=y[,2] %>% setNames(y[,1])
y2=sort(y2)
if (names(y2)[1]=="NK2A") {
  root=which(dm_coor[,1]==max(dm_coor[,1]))
} else {
  root=which(dm_coor[,1]==max(dm_coor[,1]))   
}

#dpt_d <- DPT(dm, tips=root)

pt<- dpt_d[[paste0("DPT",root)]]
pt[pt>quantile(pt,0.99,na.rm=TRUE)]=NA
pt[pt<quantile(pt,0.01,na.rm=TRUE)]=NA

df=data.frame("pseudotime"=pt,cluster=group)
rownames(df)=names(dm$DC1)  

#for colors 
pt_color<- dpt_d[[paste0("DPT",root)]]
pt_color[pt_color>quantile(pt_color,0.99,na.rm=TRUE)]=quantile(pt_color,0.99,na.rm=TRUE)
pt_color[pt_color<quantile(pt_color,0.01,na.rm=TRUE)]=quantile(pt_color,0.01,na.rm=TRUE)

clustering_d=metadata[rownames(dm_coor),"FirstClust"] %>%
  as.character() %>% setNames(rownames(dm_coor))




#Other tries
#plot(dpt_d, root= 1, col= col_firstCl[    clustering_d])
#plot(dpt_d, pch=20, col= unname(col_firstCl[    clustering_d]), frame.plot = TRUE)  + theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())
#plot(dpt_d, pch=20, col= unname(col_firstCl[    clustering_d]), frame.plot = TRUE )  + theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())

#plot(dpt_d, root= 1, col= col_firstCl[    clustering_d])
#plot(dpt_d, col= col_firstCl[    clustering_d], divide = 3, dcs = c(-1,-3,2), pch = 20)
#plot(dpt_d,  divide = 3, dcs = c(-1,-3,2), pch = 20)


# 2D plots with GGPLOT
pt<- dpt_d[[paste0("DPT",root)]]

#Reverse DC1
dm_coor[,"DC1"] = -dm_coor[,"DC1"]


object_d=subset(object,cells=rownames(dm_coor))
object_d[["dm"]]=CreateDimReducObject(dm_coor, key="DC")
object_d$pt=pt

df2=data.frame(dm_coor[,1:3],pseudotime=pt,pseudotime2=pt_color)

# DC 1 &2
p1=DimPlot(object_d,reduction="dm", group.by="FirstClust",
           cols = col_firstCl) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())  + ggtitle("")


p2=ggplot(df2, aes( DC1, DC2,col = pseudotime2)) +
  geom_point() +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_color_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))


p=p1+p2
print(p)


# Pseudotime Plot

p3 = ggplot(df, aes( pseudotime, cluster,col = cluster)) +  
  ggbeeswarm::geom_quasirandom( alpha=0.7,
                                cex=1,
                                show.legend = FALSE,
                                
                                groupOnX=FALSE)+ylab("")  + 
  ggtitle("Pseudotime in Dataset4") + theme_light()+
  stat_summary(
    aes(group = cluster), fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", color = "black", width = 0.7, lwd = 0.2,
    
    # add this bit here to your stat_summary function
    position=position_dodge(width=0.75)
  ) +
  scale_color_manual(values=col_firstCl)

p3

#Same but split pos / neg


pt<- dpt_d[[paste0("DPT",root)]]
pt[pt>quantile(pt,0.99,na.rm=TRUE)]=NA
pt[pt<quantile(pt,0.01,na.rm=TRUE)]=NA

group=factor(metadata[rownames(dm_coor),"FirstClust"],
             levels=c("NK2A","NK2B","NK1A","NK1B","NK1C","NK3"))


CMV=metadata[rownames(dm_coor),"Dataset2"]
CMV=gsub("Dataset4_","",CMV)

df=data.frame("pseudotime"=pt,cluster=group, CMV=   CMV)
df$cluster2=paste0(df$cluster,"_",df$CMV)
df$cluster2=factor(df$cluster2,
                   levels=c("NK2A_neg","NK2A_pos",
                            "NK2B_neg","NK2B_pos",
                            "NK1A_neg","NK1A_pos",
                            "NK1B_neg","NK1B_pos",
                            "NK1C_neg","NK1C_pos",
                            "NK3_neg","NK3_pos")
)
rownames(df)=rownames(dm_coor)


#PseudoTime Plot separating POS and NEG
p4 = ggplot(df, aes( pseudotime, cluster2,col = cluster)) +  
  ggbeeswarm::geom_quasirandom( alpha=0.7,
                                cex=1,
                                show.legend = FALSE,
                                #size=0.1,
                                width = 0.4,
                                groupOnX=FALSE)+ylab("")  + 
  ggtitle("Pseudotime in Dataset4") + theme_light()+
  stat_summary(
    aes(group = cluster), fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", color = "black", width = 0.7, lwd = 0.2,
    
    # add this bit here to your stat_summary function
    position=position_dodge(width=0.75)
  ) +
  scale_color_manual(values=col_firstCl) 

p4


p3

p4

#Saving Ouputs in HD

pdf("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/Figures/AllDatasets_NK1andNK2/DiffusionMap_and_DMPT.pdf",  width = 50, height = 15,  units = "cm", res=600)
p
dev.off()

svg("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/Figures/AllDatasets_NK1andNK2/DiffusionMap_and_DMPT.svg",  width = 50, height = 15,  units = "cm", res=600)
p
dev.off()


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/Figures/AllDatasets_NK1andNK2/DiffusionMap_and_DMPT.png", width = 50, height = 15,  units = "cm", res=600 )
p
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/Figures/AllDatasets_NK1andNK2/PseudotimePlot_All.png", width = 30, height = 10,  units = "cm", res=600 )
p3
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/Figures/AllDatasets_NK1andNK2/PseudotimePlot_All_splitPosNeg.png", width = 30, height = 10,  units = "cm", res=600 )
p4
dev.off()


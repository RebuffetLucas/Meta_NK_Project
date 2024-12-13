#Data description / Global BarPlots


library(SeuratDisk)

#Load Data
PBMC_Blood = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK_blood.h5seurat")
PBMC_All = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5seurat")

table(PBMC_All$meta_histology)

table(PBMC_All$meta_tissue)

table(PBMC_All$meta_histology, PBMC_All$meta_tissue)


PBMC_Blood_Bright = subset(PBMC_All, subset = Majortype == "CD56brightCD16lo" & meta_tissue ==  "Blood" )
PBMC_Blood_Dim = subset(PBMC_All, subset = Majortype == "CD56dimCD16hi" & meta_tissue ==  "Blood" )
PBMC_Blood_DP = subset(PBMC_All, subset = Majortype == "CD56brightCD16hi" & meta_tissue ==  "Blood" )

PBMC_Blood_All = subset(PBMC_All, subset = meta_tissue ==  "Blood" )

levels(PBMC_Blood$celltype) =  levels(PBMC_All$celltype)
  
PBMC = PBMC_Blood

#Blood
table(PBMC_Blood$celltype, PBMC_Blood$meta_histology)
table(PBMC_Blood$meta_histology)

table(PBMC_Blood_All$celltype)
table(PBMC_Blood_All$meta_histology)


DimPlot(PBMC_Blood, group.by = "celltype")
DimPlot(PBMC_Blood, group.by = "Majortype")

#All
table(PBMC_Blood_All$meta_histology)
table(PBMC_Blood_All$meta_tissue)
table(PBMC_Blood_All$celltype)
table(PBMC_Blood_All$meta_histology)
DimPlot(PBMC_Blood_All, group.by = "celltype", split.by = "Majortype")
DimPlot(PBMC_Blood_All, group.by = "Majortype")



#Pie Chart proportions
#Get a pie chart of the proportions of each cluster in the pop

#For Bright
df= as.data.frame(table(PBMC$Majortype))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

p5 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    

p6 = DimPlot(PBMC_Blood, group.by = "Majortype")

p5 + p6 


#For Dim

PBMC_Blood_Bright$celltype = droplevels(PBMC_Blood_Bright$celltype)
df= as.data.frame(table(PBMC_Blood_Bright$celltype))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

df=df[df$Freq>0.01,]

p7 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    



PBMC_Blood_Dim$celltype = droplevels(PBMC_Blood_Dim$celltype)
df= as.data.frame(table(PBMC_Blood_Dim$celltype))
colnames(df)[1] = "Cluster"
colnames(df)[2] = "Number"
df$Freq = df$Number / sum(df$Number)

df=df[df$Freq>0.01,]

p8 = ggplot(df, aes(x="", y= Freq , fill= Cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() +   geom_text_repel(aes(label = paste0(round(Freq*100), "%") ),   position = position_stack(vjust = 0.5), size = 6) + theme(text = element_text(size = 20))    

p7 + p8


#BarPlot Blood Composition / cancer type



p9 = ggplot(PBMC_Blood_Dim@meta.data, aes(x=meta_histology, fill= celltype)) + geom_bar(position="fill", width = 0.8) + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))
p10 = ggplot(PBMC_Blood_Bright@meta.data, aes(x=meta_histology, fill= celltype)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  
p11 = ggplot(PBMC_Blood_All@meta.data, aes(x=meta_histology, fill= Majortype)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  

p9 + p10

p11

p12 = DimPlot(PBMC_All , group.by = "celltype", split.by = "Majortype")
p13 = DimPlot(PBMC_All , group.by = "meta_tissue", split.by = "Majortype")

p14 = DimPlot(PBMC_Blood_Dim, group.by = "celltype")
p15 = DimPlot(PBMC_Blood_Bright, group.by = "celltype")

p12
p13 
p14 + p13

p15 = DimPlot(PBMC_All , group.by = "Majortype")

p15

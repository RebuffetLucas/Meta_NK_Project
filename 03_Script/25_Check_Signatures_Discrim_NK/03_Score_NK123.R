#Applying scores NK1 , NK2 and NK3:

#Extract signatures
Markers_NK_123_CITEseq = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds")

Markers_NK_123_CITEseq %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All

top_All %>% filter(cluster== "NK_1") -> Top_NK1
top_All %>% filter(cluster== "NK_2") -> Top_NK2
top_All %>% filter(cluster== "NK_3") -> Top_NK3

#have a look at the gene list:
list(Top_NK1$gene)
list(Top_NK2$gene)
list(Top_NK3$gene)

#Apply scores
PBMC = AddModuleScore(PBMC, features = list(Top_NK1$gene), pool= NULL , name= "NK1score" , seed=19)
PBMC = AddModuleScore(PBMC, features = list(Top_NK2$gene), pool= NULL , name= "NK2score" , seed=19)
PBMC = AddModuleScore(PBMC, features = list(Top_NK3$gene), pool= NULL , name= "NK3score" , seed=19)

#Convert to factors:
PBMC$celltype.l1 = as.factor(PBMC$celltype.l1)
PBMC$celltype.l2 = as.factor(PBMC$celltype.l2)
PBMC$celltype.l3 = as.factor(PBMC$celltype.l3)




#Dataviz and save

#Level 1:
PBMC = SetIdent(PBMC,  value= "celltype.l1" )

#NK1
p1 = VlnPlot(PBMC,features = "NK1score1", pt.size = 0, sort = "decreasing") & ggtitle("NK1 score") & stat_summary(fun.data=data_summary,color="black")
  print(p1)
  
png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK1score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


#NK2  
p1 = VlnPlot(PBMC,features = "NK2score1", pt.size = 0, sort = "decreasing") & ggtitle("NK2 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK2score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


#NK3
p1 = VlnPlot(PBMC,features = "NK3score1", pt.size = 0, sort = "decreasing") & ggtitle("NK3 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK3score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()



#Level 2:
PBMC = SetIdent(PBMC,  value= "celltype.l2" )

#NK1
p1 = VlnPlot(PBMC,features = "NK1score1", pt.size = 0, sort = "decreasing") & ggtitle("NK1 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK1score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


#NK2  
p1 = VlnPlot(PBMC,features = "NK2score1", pt.size = 0, sort = "decreasing") & ggtitle("NK2 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK2score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


#NK3
p1 = VlnPlot(PBMC,features = "NK3score1", pt.size = 0, sort = "decreasing") & ggtitle("NK3 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK3score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()



#Level 3:
PBMC = SetIdent(PBMC,  value= "celltype.l3" )

#NK1
p1 = VlnPlot(PBMC,features = "NK1score1", pt.size = 0, sort = "decreasing") & ggtitle("NK1 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK1score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


#NK2  
p1 = VlnPlot(PBMC,features = "NK2score1", pt.size = 0, sort = "decreasing") & ggtitle("NK2 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK2score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()


#NK3
p1 = VlnPlot(PBMC,features = "NK3score1", pt.size = 0, sort = "decreasing") & ggtitle("NK3 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_NK3score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()




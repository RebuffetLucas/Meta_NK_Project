#Applying scores NK1 , NK2 and NK3:
#6 pops
PBMC= reference


#Extract the genes
Markers_6_clusters = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/All_Markers_6clusters.rds")


Markers_6_clusters %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All_6cl



top_All_6cl %>% filter(cluster== "NK1A") -> Top_NK1A_6cl
top_All_6cl %>% filter(cluster== "NK1B") -> Top_NK1B_6cl
top_All_6cl %>% filter(cluster== "NK1C") -> Top_NK1C_6cl
top_All_6cl %>% filter(cluster== "NKint") -> Top_NKint_6cl
top_All_6cl %>% filter(cluster== "NK2") -> Top_NK2_6cl
top_All_6cl %>% filter(cluster== "NK3") -> Top_NK3_6cl

#Liste des gènes
list(Top_NK1A_6cl$gene) #Carefull only 10 genes here!
list(Top_NK1B_6cl$gene)
list(Top_NK1C_6cl$gene)
list(Top_NKint_6cl$gene)
list(Top_NK2_6cl$gene)
list(Top_NK3_6cl$gene)


#Apply scores
PBMC = AddModuleScore(PBMC, features = list(Top_NK1A_6cl$gene), pool= NULL , name= "NK1Ascore" , seed=19)
PBMC = AddModuleScore(PBMC, features = list(Top_NK1B_6cl$gene), pool= NULL , name= "NK1Bscore" , seed=19)
PBMC = AddModuleScore(PBMC, features = list(Top_NK1C_6cl$gene), pool= NULL , name= "NK1Cscore" , seed=19)

PBMC = AddModuleScore(PBMC, features = list(Top_NKint_6cl$gene), pool= NULL , name= "NKintscore" , seed=19)
PBMC = AddModuleScore(PBMC, features = list(Top_NK2_6cl$gene), pool= NULL , name= "NK2score" , seed=19)
PBMC = AddModuleScore(PBMC, features = list(Top_NK3_6cl$gene), pool= NULL , name= "NK3score" , seed=19)




#Dataviz and figure saving

#Level 1:
PBMC = SetIdent(PBMC,  value= "celltype.l1" )

#NK1A
p1 = VlnPlot(PBMC,features = "NK1Ascore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1A score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1A","score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK1B
p1 = VlnPlot(PBMC,features = "NK1Bscore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1B score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1B","score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK1C
p1 = VlnPlot(PBMC,features = "NK1Cscore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1C score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1C","score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NKint
p1 = VlnPlot(PBMC,features = "NKintscore1", pt.size = 0, sort = "decreasing") & ggtitle("NKint score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NKint","score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK2
p1 = VlnPlot(PBMC,features = "NK2score1", pt.size = 0, sort = "decreasing") & ggtitle("NK2 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK2","score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK3
p1 = VlnPlot(PBMC,features = "NK3score1", pt.size = 0, sort = "decreasing") & ggtitle("NK3 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK3","score_OnCD45pos_", "Level1" ,"_incr.png"), width = 25, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()




#Level 2:
PBMC = SetIdent(PBMC,  value= "celltype.l2" )

#NK1A
p1 = VlnPlot(PBMC,features = "NK1Ascore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1A score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1A","score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK1B
p1 = VlnPlot(PBMC,features = "NK1Bscore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1B score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1B","score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK1C
p1 = VlnPlot(PBMC,features = "NK1Cscore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1C score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1C","score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NKint
p1 = VlnPlot(PBMC,features = "NKintscore1", pt.size = 0, sort = "decreasing") & ggtitle("NKint score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NKint","score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK2
p1 = VlnPlot(PBMC,features = "NK2score1", pt.size = 0, sort = "decreasing") & ggtitle("NK2 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK2","score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK3
p1 = VlnPlot(PBMC,features = "NK3score1", pt.size = 0, sort = "decreasing") & ggtitle("NK3 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK3","score_OnCD45pos_", "Level2" ,"_incr.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()










#Level 3:
PBMC = SetIdent(PBMC,  value= "celltype.l3" )

#NK1A
p1 = VlnPlot(PBMC,features = "NK1Ascore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1A score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1A","score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK1B
p1 = VlnPlot(PBMC,features = "NK1Bscore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1B score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1B","score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK1C
p1 = VlnPlot(PBMC,features = "NK1Cscore1", pt.size = 0, sort = "decreasing") & ggtitle("NK1C score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK1C","score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NKint
p1 = VlnPlot(PBMC,features = "NKintscore1", pt.size = 0, sort = "decreasing") & ggtitle("NKint score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NKint","score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK2
p1 = VlnPlot(PBMC,features = "NK2score1", pt.size = 0, sort = "decreasing") & ggtitle("NK2 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK2","score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()

#NK3
p1 = VlnPlot(PBMC,features = "NK3score1", pt.size = 0, sort = "decreasing") & ggtitle("NK3 score") & stat_summary(fun.data=data_summary,color="black")
print(p1)

png(file=paste0(OUPUT_FIG_FILE,"/VlnPlot_","NK3","score_OnCD45pos_", "Level3" ,"_incr.png"), width = 60, height = 15,  units = "cm", res=600 )
print(p1)
dev.off()






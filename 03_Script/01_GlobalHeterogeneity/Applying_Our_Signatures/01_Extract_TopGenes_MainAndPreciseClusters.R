#Extracting our signatures For Dim/Bright/Adapt and for NK1A,NK1B, ...

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

PBMC$MainClust = PBMC$FirstClust

levels(PBMC$MainClust) = c("NK1" , "NK3" , "NK1", "NK1" , "NK2" , "NK2")

PBMC=SetIdent(PBMC, value = "MainClust")

DimPlot(PBMC)

#Look at markers previously calculated
#Find Markers
All_Markers = FindAllMarkers(PBMC , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


#Look at markers like Constance
markers_overCl=All_Markers

#Make a black_list
mito.genes <- grep(pattern = "^MT-", x = rownames(PBMC), value = TRUE)
ribo.genes <- grep(pattern = "^RPS", x = rownames(PBMC), value = TRUE)
ribo.genes2 <- grep(pattern = "^RPL", x = rownames(PBMC), value = TRUE)

black_list= Reduce(union, list(mito.genes, ribo.genes, ribo.genes2))

#Sort genes of interest
All_Markers %>% 
  group_by(cluster) %>%
  filter(gene %in% rownames(PBMC) ) %>%
  filter(!gene %in% black_list ) %>%
  filter(p_val_adj<0.05 & avg_log2FC>0 ) %>%
  top_n(n=FINDMARKERS_SHOWTOP , wt = avg_log2FC) -> top10


table(top10$cluster)


List_Top10=list()

for (i in levels(top10$cluster)){
  top10 %>% 
    filter(cluster== i ) -> top10_cluster
  print(top10_cluster$gene)
  List_Top10  = append(List_Top10, data.frame(top10_cluster$gene))
}

names(List_Top10) = levels(top10$cluster)


List_Top10= lapply(List_Top10, as.data.frame)

#Change colnames
for (i in names(List_Top10)){
  print(i)
  colnames(List_Top10[[i]]) = i
}

saveRDS(List_Top10, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds" )


#For Bright:
PBMC2=PBMC
PBMC= subset(PBMC2, idents="NK2" )
PBMC= SetIdent(PBMC, value = "FirstClust")

#rerun the script with PBMC= Bright




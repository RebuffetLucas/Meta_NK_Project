
A = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg/de_up_adaptive.csv")
B = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg/de_up_brights.csv")
C = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg/de_up_cd57_pos.csv")
D = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg/de_up_kir_pos.csv")
E = read.csv("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg/de_up_nkg2a_pos.csv")

PBMC = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

List_Malmberg=list(A$X , B$X, C$X,D$X,E$X)
names(List_Malmberg)= c("Adaptive","Bright","CD57pos","KIRpos","NKG2Apos")

List_df_Malmberg= lapply(List_Malmberg, as.data.frame)


#Change colnames
for (i in names(List_df_Malmberg)){
  print(i)
  colnames(List_df_Malmberg[[i]]) = i
}


saveRDS(List_df_Malmberg, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg.rds" )


#Only names in RowNames AND same length (here TOP10) to be able to compare and attribute maximum
List_df_Malmberg= lapply(List_df_Malmberg, FUN=function(x) intersect(rownames(PBMC), as.character(unlist(x))))
min_genes= min(unlist(lapply(List_df_Malmberg, FUN=function(x) length(x))))
List_df_Malmberg= lapply(List_df_Malmberg, FUN=function(x) head(x , n= min_genes))
List_df_Malmberg= lapply(List_df_Malmberg, as.data.frame)
#Change colnames
for (i in names(List_df_Malmberg)){
  print(i)
  colnames(List_df_Malmberg[[i]]) = i
}


saveRDS(List_df_Malmberg, "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg_TOP10.rds" )




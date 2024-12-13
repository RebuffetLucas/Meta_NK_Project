List_query = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/22_Comparing_Colonna_paper/Colonna_Transfer_Label_SCT/data.query_all_cells_Healthy_Tissues_Colonna_SCTnorm_AllMAppedatONCE.rds")

List_query = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/22_Comparing_Colonna_paper/Colonna_Transfer_Label_SCT/List_query_IterativelyGenerated_SCTNorm.rds")




QUERY_DATASETS = names(List_query)

data.query = merge(List_query[QUERY_DATASETS][[1]] , y = List_query[QUERY_DATASETS][-1]) #Merge all the query together


table(data.query$predicted.id , data.query$meta_histology)

table(data.query$predicted.id , data.query$SecondClust )


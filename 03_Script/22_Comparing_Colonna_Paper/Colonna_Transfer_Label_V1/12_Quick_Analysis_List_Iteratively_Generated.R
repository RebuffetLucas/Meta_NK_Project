List_query = readRDS(file=paste0(PATH_ANALYSIS_OUTPUT,"/List_query_ItreativelyGenerated.rds"))


QUERY_DATASETS = names(List_query)

data.query = merge(List_query[QUERY_DATASETS][[1]] , y = List_query[QUERY_DATASETS][-1]) #Merge all the query together


table(data.query$predicted.id , data.query$meta_histology)

table(data.query$predicted.id , data.query$SecondClust )


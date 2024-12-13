#############################################################
## Find descriptions of interest####
###############################################################

library(DOSE)
TERMS_TO_REMOVE= "neuron|virus|symbio|host|folding|folded|topologically|viral|bacterial|lipopolysaccharide|RNA|ribosome|ribo|translation|transcription|ubiquitin|amide|chromosome|nucleosom|nucleoplasm|conjugation|mitochon|telomere|chromatin|organelle|erythrocyte|number of cells|B cell|arsenic|Cajal|proteolysis"

TERMS_TO_REMOVE= NULL



## @knitr clusterComp_Visualisation


for(j in 1:9){
  cat("\n\n") # Required for '.tabset'
  cat("#### ",names(CompareCLUST[j]), "\n\n ") # tabset title
  
  enrichplot::dotplot(CompareCLUST[[j]], font.size = 10, label_format = 90, includeAll=FALSE) + 
    theme(axis.text.x=element_text(size = 10), axis.text.y=element_text(size = 10), 
          legend.text =  element_text(size = 10), legend.title =  element_text(size = 10))
  
  cat(" \n \n") # Required for '.tabset'
}


#Print Only One of them:
COMP_OF_INTEREST = 3
enrichplot::dotplot(CompareCLUST[[COMP_OF_INTEREST]], font.size = 10, label_format = 90, showCategory=30 ,includeAll=TRUE) + 
  theme(axis.text.x=element_text(size = 10), axis.text.y=element_text(size = 10), 
        legend.text =  element_text(size = 10), legend.title =  element_text(size = 10))



GO_results = CompareCLUST$GO

Go_results_CompareClust = GO_results@compareClusterResult

#Have a look at Shared Descriptions

extract_duplicates <- function(lst) {
  table_lst <- table(lst)
  duplicates <- names(table_lst[table_lst > 1])
  return(duplicates)
}


Shared_Descriptions = extract_duplicates(Go_results_CompareClust$Description)
GO_to_plot



#Remove unwanted terms
all_Descriptions = unique(Go_results_CompareClust$Description)

black_List2 <- data.frame(all_Descriptions) %>%
  filter( grepl(TERMS_TO_REMOVE, all_Descriptions))

black_List2= black_List2$all_Descriptions

black_List = black_List2


Go_results_CompareClust2 = filter(Go_results_CompareClust,  !(Description %in% black_List))


Go_results_CompareClust2$GeneRatio = parse_ratio(Go_results_CompareClust2$GeneRatio)

Go_results_CompareClust2 %>%
  group_by(Cluster) %>%
  top_n(n=40, wt = -p.adjust) -> top10

#Extract in Top10 the most interesting 
unique(top10$Description)[100:110]





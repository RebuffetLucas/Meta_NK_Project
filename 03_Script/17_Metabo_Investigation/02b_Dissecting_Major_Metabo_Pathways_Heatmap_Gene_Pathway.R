#Create A umap OF THE USUAL SUSPECTS

library(rlist)
library(reshape2)


#Get the list of genes present in at list X percent of the cells
MIN_PCT_GENE= 10

MIN_ACCEPTABLE_GENES = 2 # We will kepp pathway with number of genes  >= 


PBMC= SetIdent(PBMC, value = "FirstClust")

#List of list of Genes for each Pathway

List_Markers = unique(Genes_To_Test)
List_of_Lists= Genes_To_Test2
names(List_of_Lists) = Summary_Pathways$Annot




#Genes that appear in at list  MIN_PCT_GENE% of the cells
perc_exp <- DotPlot(PBMC, features=List_Markers, group.by="FirstClust")$data[, c("features.plot", "id", "pct.exp")]

perc_exp %>%
  filter(pct.exp > MIN_PCT_GENE ) -> Acceptable_Genes

Acceptable_Genes = unique(as.character(Acceptable_Genes$features.plot))

  #Quick Check

length(List_Markers)
length(Acceptable_Genes)

#Keep only these acceptable genes
List_total= intersect(List_Markers , Acceptable_Genes  )

List_of_Lists=  lapply(List_of_Lists , function(x) intersect(x, Acceptable_Genes))

str(Genes_To_Test2)
str(List_of_Lists)

#Keep only the list with at least one gene
List_of_Lists2 = Filter(function(x) length(x) >= MIN_ACCEPTABLE_GENES, List_of_Lists)



#Calculate Average expression within each cluster
PBMC$FirstClust= droplevels(PBMC$FirstClust)
PBMC= SetIdent(PBMC, value = "FirstClust")


#Try with ComplexHeatmap:
cluster.averages <- AverageExpression(PBMC, group.by = "FirstClust", features = List_total , slot= "data")
mat = cluster.averages$RNA
#mat.scaled= t(scale(t(mat)))
#mat.scaled = (mat-apply(mat, 1, mean) / apply(mat, 1, sd))
mat.scaled= t(apply(mat, MARGIN = 1, FUN = function(X) (X - mean(X))/sd(X)))
col_fun = circlize::colorRamp2(c(-2,-1, 0, 1,2), rev(brewer.pal(n=5,name="RdBu")))
#col_fun = circlize::colorRamp2(c(-4,-2, -1, -0.5, 0, 0.5, 1, 2, 4), rev(brewer.pal(n=9,name="RdBu")))
#col_fun = circlize::colorRamp2(c(-4,-3 , -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4), rev(brewer.pal(n=11,name="RdBu")))

SIZE_ROW = 8

setdiff( names(List_of_Lists) ,  names(List_of_Lists2)) #Pathways that did not pass the filter




#p2 = ComplexHeatmap::Heatmap(mat.scaled[intersect(List_of_Lists2[[2]], rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col= col_fun , border = TRUE, show_heatmap_legend= FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))

  #Make groups of heatmaps that are more or less equal in terms of size:
splitList <- function(listOfLists) {
  # Initialize variables
  groups <- list()
  currentGroup <- list()
  currentLength <- 0
  
  # Iterate over the list of lists
  for (name in names(listOfLists)) {
    sublist <- listOfLists[[name]]
    sublistLength <- length(sublist)
    
    # Check if adding the current sublist exceeds the length limit
    if ((currentLength + sublistLength) > 100) {
      # Save the current group and start a new one
      groups[[length(groups) + 1]] <- currentGroup
      currentGroup <- list()
      currentLength <- 0
    }
    
    # Add the sublist to the current group
    currentGroup[[name]] <- sublist
    currentLength <- currentLength + sublistLength
  }
  
  # Add the last group if it's not empty
  if (length(currentGroup) > 0) {
    groups[[length(groups) + 1]] <- currentGroup
  }
  
  return(groups)
}


List_of_Lists3 = splitList(List_of_Lists2)

#Create one last with Glycolysis + TCA + Cystein 
#Add one with TCA + Glycolysis + Cystein
# Extracting the specific lists
glycolysis_gluconeogenesis <- List_of_Lists3[[2]][[1]]
citrate_cycle <- List_of_Lists3[[1]][[7]]
cysteine_methionine_metabolism <- List_of_Lists3[[1]][[8]]

# Creating the new sixth sublist with proper naming
List_of_Lists3[[6]] <- list(
  `Glycolysis / Gluconeogenesis` = glycolysis_gluconeogenesis,
  `Citrate cycle (TCA cycle)` = citrate_cycle,
  `Cysteine and methionine metabolism` = cysteine_methionine_metabolism
)

List_of_Lists3



# Initialize a list to store heatmap objects
heatmap_list <- list()

for (i in 1:length(List_of_Lists3)) {
  # Initialize an empty plot list for each i
  plot_list <- list()
  
  for (j in 1:length(List_of_Lists3[[i]])) {
    pathway_genes = names(List_of_Lists3[[i]][[j]])
    gene_set = List_of_Lists3[[i]][[j]]
    
    # Check if the gene set is not empty and has matching rownames in cluster.averages$RNA
    if (length(gene_set) > 0 && any(gene_set %in% rownames(cluster.averages$RNA))) {
      # Create a heatmap with or without the legend, depending on its position
      if (j == 1) {
        # Include legend for the first heatmap
        heatmap <- ComplexHeatmap::Heatmap(mat.scaled[intersect(gene_set, rownames(cluster.averages$RNA)),], 
                                           cluster_columns = FALSE, 
                                           col = col_fun, 
                                           border = TRUE, 
                                           row_names_gp = gpar(fontsize = SIZE_ROW),
                                           heatmap_legend_param = list(title = "Z score", 
                                                                       direction = "horizontal", 
                                                                       legend_width = unit(10, "cm")))
      } else {
        # Exclude legend for subsequent heatmaps
        heatmap <- ComplexHeatmap::Heatmap(mat.scaled[intersect(gene_set, rownames(cluster.averages$RNA)),], 
                                           cluster_columns = FALSE, 
                                           col = col_fun, 
                                           border = TRUE, 
                                           row_names_gp = gpar(fontsize = SIZE_ROW),
                                           show_heatmap_legend = FALSE)
      }
      
      # Add the heatmap to the plot list
      plot_list[[j]] <- heatmap
    }
  }
  
  # Combine all heatmaps for the current i into one and store in the main list
  if (length(plot_list) > 0) {
    combined_heatmap <- Reduce(`%v%`, plot_list)
    heatmap_list[[i]] <- combined_heatmap
  }
}


p1 =draw(heatmap_list[[1]], heatmap_legend_side = "top")

p2 = draw(heatmap_list[[2]], heatmap_legend_side = "top")

p3 = draw(heatmap_list[[3]], heatmap_legend_side = "top")

p4 =draw(heatmap_list[[4]], heatmap_legend_side = "top")

p5= draw(heatmap_list[[5]], heatmap_legend_side = "top")

p6 = draw(heatmap_list[[6]], heatmap_legend_side = "top") #For TCA + Glycolysis + Cystein












#Save Them


png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Heatmap_Genes_Pathways/V2/HeatMap_p1.png", width = 18, height = 30,  units = "cm", res=600 )
p1
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Heatmap_Genes_Pathways/V2/HeatMap_p2.png", width = 18, height = 30,  units = "cm", res=600 )
p2
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Heatmap_Genes_Pathways/V2/HeatMap_p3.png", width = 18, height = 30,  units = "cm", res=600 )
p3
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Heatmap_Genes_Pathways/V2/HeatMap_p4.png", width = 18, height = 30,  units = "cm", res=600 )
p4
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Heatmap_Genes_Pathways/V2/HeatMap_p5.png", width = 18, height = 30,  units = "cm", res=600 )
p5
dev.off()

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/17_Metabo_Investigation/Deeper_Investigation/Heatmap_Genes_Pathways/V2/HeatMap_p6_TCA_Glyco_Cysteine.png", width = 18, height = 30,  units = "cm", res=600 )
p6
dev.off()



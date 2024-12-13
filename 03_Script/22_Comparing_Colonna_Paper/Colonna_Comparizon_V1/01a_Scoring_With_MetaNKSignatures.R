## @knitr Scoring_With_MetaNKSignatures

#def function

  # +- std
data_summary <- function(x) {
   m <- median(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}



#Importing the rds files of Colonna

cat(" \n \n")
cat("# Comparizon with Colonna's paper")
cat(" \n \n")



cat(" \n \n")
cat("## Scoring our data with ILC and ILCP signatures (from their Fig 1a)")
cat(" \n \n")

PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")
PBMC_Meta_NK_V2 = SetIdent(PBMC_Meta_NK_V2 , value = "FirstClust")

#Scoring with our signatures of 3 main pops
Markers_Fig1a = read.csv("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/markers_Fig1a.csv")


Markers_Fig1a %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_Fig1a

top_Fig1a %>% filter(cluster== "ILC?") -> Top_ILC
top_Fig1a %>% filter(cluster== "ILCP") -> Top_ILCP
top_Fig1a %>% filter(cluster== "ILC2") -> Top_ILC2


PBMC_Meta_NK_V2 = AddModuleScore( PBMC_Meta_NK_V2 , features = list(Top_ILC$gene), pool= NULL ,name= "ILC_notsure" , seed=19) 
PBMC_Meta_NK_V2 = AddModuleScore( PBMC_Meta_NK_V2 , features = list(Top_ILCP$gene) , pool= NULL ,name= "ILCP" , seed=19) 
PBMC_Meta_NK_V2 = AddModuleScore( PBMC_Meta_NK_V2 , features = list(Top_ILC2$gene) , pool= NULL ,name= "ILC2" , seed=19) 



cat(" \n \n")
p1 =DimPlot(PBMC_Meta_NK_V2, pt.size = PT_SIZE, cols = palette) & ggtitle("Meta NK")
print(p1)
cat(" \n \n")

cat(" \n \n")
p2= FeaturePlot(PBMC_Meta_NK_V2, features="ILC_notsure1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("ILC_notsure")   & NoAxes()
p3= FeaturePlot(PBMC_Meta_NK_V2, features="ILCP1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("ILCP")   & NoAxes()

cat(" \n \n")
print(p2)
cat(" \n \n")

cat(" \n \n")
print(p3)
cat(" \n \n")


p4= VlnPlot(PBMC_Meta_NK_V2, features="ILC_notsure1" , pt.size = PT_SIZE_VLN, y.max = 1.2, cols = palette) & ggtitle("ILC") & stat_summary(fun.data=data_summary,color="black")

p5= VlnPlot(PBMC_Meta_NK_V2, features="ILCP1" , pt.size = PT_SIZE_VLN, y.max=1.2, cols = palette) & ggtitle("ILCP") & stat_summary(fun.data=data_summary,color="black")

p6= VlnPlot(PBMC_Meta_NK_V2, features="ILC21" , pt.size = PT_SIZE_VLN, y.max=1., cols = palette) & ggtitle("ILC2") & stat_summary(fun.data=data_summary,color="black")



cat(" \n \n")
print(p4)
cat(" \n \n")

cat(" \n \n")
print(p5)
cat(" \n \n")

cat(" \n \n")
print(p6)
cat(" \n \n")


FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/VlnPlot_Scoring_ILC_notsure.png")

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p4)
dev.off()





FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/VlnPlot_Scoring_ILCP.png")
  
png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p5)
dev.off()

FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/VlnPlot_Scoring_ILC2.png")

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p6)
dev.off()

p6= RidgePlot(PBMC_Meta_NK_V2, features="ILC_notsure1" ) & ggtitle("ILC_notsure")
p7= RidgePlot(PBMC_Meta_NK_V2, features="ILCP1" ) & ggtitle("ILCP")


cat(" \n \n")
print(p6)
cat(" \n \n")

cat(" \n \n")
print(p7)
cat(" \n \n")



cat(" \n \n")
cat(" \n \n")



cat(" \n \n")
cat("## Scoring Colonna's data {.tabset .tabset-fade} \n\n")
cat(" \n \n")





###Loading the data
Co_ILCs_Blood_1a = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Blood_Fig1a.rds")
#DimPlot(Co_ILCs_Blood_1a, group.by = "final_nomenclature" )
Co_ILCs_Blood_1a = SetIdent(Co_ILCs_Blood_1a, value = "final_nomenclature")


Co_NK_Blood_1f = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Blood_Fig1f.rds")
#DimPlot(Co_NK_Blood_1f, group.by = "final_nomenclature")
Co_NK_Blood_1f = SetIdent(Co_NK_Blood_1f, value = "final_nomenclature")

Co_NK_Tonsil_2d = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Tonsil_Fig2d.rds")
#DimPlot(Co_NK_Tonsil_2d, group.by = "final_nomenclature")
Co_NK_Tonsil_2d = SetIdent(Co_NK_Tonsil_2d, value = "final_nomenclature")

Co_NK_Tonsil_2d_NK_Only = subset(Co_NK_Tonsil_2d , idents= c("8 -ILC3"), invert=TRUE )
Co_NK_Tonsil_2d_NK_Only = SCTransform(Co_NK_Tonsil_2d_NK_Only, verbose = FALSE)
Co_NK_Tonsil_2d_NK_Only$final_nomenclature = droplevels(Co_NK_Tonsil_2d_NK_Only$final_nomenclature)


Co_NK_Lung_3a = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Lung_Fig3a.rds")
#DimPlot(Co_NK_Lung_3a , group.by = "final_nomenclature")
Co_NK_Lung_3a$final_nomenclature = as.factor(Co_NK_Lung_3a$final_nomenclature)
Co_NK_Lung_3a = SetIdent(Co_NK_Lung_3a, value = "final_nomenclature")

Co_NK_Lung_3a_NK_Only = subset(Co_NK_Lung_3a, idents = c(  "ILC2_ILC3" ), invert = TRUE)
#Co_NK_Lung_3a_NK_Only = subset(Co_NK_Lung_3a, idents = c("ILC1",  "ILC2_ILC3" ), invert = TRUE)
Co_NK_Lung_3a_NK_Only=SCTransform(Co_NK_Lung_3a_NK_Only,verbose = FALSE)
Co_NK_Lung_3a_NK_Only$final_nomenclature =droplevels(Co_NK_Lung_3a_NK_Only$final_nomenclature )




Co_ILCs_Intestines_4k = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/IEL_LPL_Fig4k.rds")
#DimPlot(Co_ILCs_Intestines_4k , group.by = "final_nomenclature")
Co_ILCs_Intestines_4k = SetIdent(Co_ILCs_Intestines_4k, value = "final_nomenclature")


levels(Co_ILCs_Intestines_4k$final_nomenclature)
Co_ILCs_Intestines_4k_NK_Only = subset(Co_ILCs_Intestines_4k , idents = c("ILC3"  ,   "ILC3p",     "ILC3-ILC1" ) , invert=TRUE)
#Co_ILCs_Intestines_4k_NK_Only = subset(Co_ILCs_Intestines_4k , idents = c("ILC3"  ,  "ILC1" ,  "ILC3p",     "ILC3-ILC1" ) , invert=TRUE)
Co_ILCs_Intestines_4k_NK_Only= SCTransform(Co_ILCs_Intestines_4k_NK_Only,verbose = FALSE)
Co_ILCs_Intestines_4k_NK_Only$final_nomenclature = droplevels( Co_ILCs_Intestines_4k_NK_Only$final_nomenclature)
#DimPlot(Co_ILCs_Intestines_4k_NK_Only)


Co_ILCs_All_Tissues_5a = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/AllTissues_Fig5a_b_c.rds")
Co_ILCs_All_Tissues_5a = SCTransform(Co_ILCs_All_Tissues_5a, verbose = FALSE) #SCTransform need to be rerun for this one
Co_ILCs_All_Tissues_5a$tissue_chief = as.factor(Co_ILCs_All_Tissues_5a$tissue_chief)
Co_ILCs_All_Tissues_5a = SetIdent(Co_ILCs_All_Tissues_5a, value = "tissue_chief")



Co_ILCs_All_Tissues_5a_NK_Only = subset( Co_ILCs_All_Tissues_5a, idents= c( "pbmc_ilcp" , "tonsil_ilcp"), invert= TRUE)
Co_ILCs_All_Tissues_5a_NK_Only= SCTransform(Co_ILCs_All_Tissues_5a_NK_Only,verbose = FALSE)
Co_ILCs_All_Tissues_5a_NK_Only$tissue_chief = droplevels(Co_ILCs_All_Tissues_5a_NK_Only$tissue_chief)

#DimPlot(Co_ILCs_All_Tissues_5a, group.by = "chief")
#DimPlot(Co_ILCs_All_Tissues_5a, group.by = "tissue")
#DimPlot(Co_ILCs_All_Tissues_5a, group.by = "tissue_chief")



#Puting them all in a gene list
Co_List_Objects = list(Co_ILCs_Blood_1a , Co_NK_Blood_1f , Co_NK_Tonsil_2d, Co_NK_Tonsil_2d_NK_Only , Co_NK_Lung_3a, Co_NK_Lung_3a_NK_Only , Co_ILCs_Intestines_4k, Co_ILCs_Intestines_4k_NK_Only ,  Co_ILCs_All_Tissues_5a, Co_ILCs_All_Tissues_5a_NK_Only)
names(Co_List_Objects) = c("Co_ILCs_Blood_1a" , "Co_NK_Blood_1f" , "Co_NK_Tonsil_2d", "Co_NK_Tonsil_2d_NK_Only", "Co_NK_Lung_3a", "Co_NK_Lung_3a_NK_Only" ,"Co_ILCs_Intestines_4k", "Co_ILCs_Intestines_4k_NK_Only", "Co_ILCs_All_Tissues_5a", "Co_ILCs_All_Tissues_5a_NK_Only")




#Main pops
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


Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK1$gene), pool= NULL ,name= "NK1" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK2$gene), pool= NULL ,name= "NK2" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK3$gene), pool= NULL ,name= "NK3" , seed=19)

Markers_NK_123_CITEseq %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_GENES_JACKARD_3pops, wt = avg_log2FC) -> top_All

#6 pops


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


Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK1A_6cl$gene), pool= NULL ,name= "NK1A_6cl" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK1B_6cl$gene), pool= NULL ,name= "NK1B_6cl" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK1C_6cl$gene), pool= NULL ,name= "NK1C_6cl" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NKint_6cl$gene), pool= NULL ,name= "NKint_6cl" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK2_6cl$gene), pool= NULL ,name= "NK2_6cl" , seed=19)
Co_List_Objects = lapply(Co_List_Objects , AddModuleScore,  features = list(Top_NK3_6cl$gene), pool= NULL ,name= "NK3_6cl" , seed=19)


Markers_6_clusters %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_GENES_JACKARD_6pops, wt = avg_log2FC) -> top_All_6cl



compteur_tour=0
for (object in Co_List_Objects){
  compteur_tour=compteur_tour+1
  names(Co_List_Objects)[compteur_tour]
  cat(" \n \n")
  cat("### ",  names(Co_List_Objects)[compteur_tour], "\n")
  cat(" \n \n")

  palette_i = hue_pal()(length(levels(object@active.ident)))
  names(palette_i) = levels(object@active.ident)
  
  cat(" \n \n")
  cat("#### Spontaneous Signatures")
  cat(" \n \n")
  
  
  p1 =DimPlot(object, pt.size = PT_SIZE, cols=palette_i, label = TRUE) & ggtitle("Clust")
  
  
  
  All_Markers = FindAllMarkers(object , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)
  
  All_Markers %>%
    group_by(cluster) %>%
    filter(!grepl("RPL|RPS|MT-", gene)) %>%
    top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top10
  
  Dot1 = DotPlot(object, features = unique(top10$gene), cols = "RdBu") + theme(axis.text.x = element_text(angle = 90))
  
  All_Markers %>%
    group_by(cluster) %>%
    filter(!grepl("RPL|RPS|MT-", gene)) %>%
    top_n(n = NUMBER_GENES_JACKARD, wt = avg_log2FC) -> top10
  
  
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  cat(" \n \n")
  print(Dot1)
  cat(" \n \n")
  
  
  cat(" \n \n")
  cat("#### Overlap of Signatures")
  cat(" \n \n")
  
  
  #Jackard and overlap Index for 6 clusters
  # Split the data frame into a list of data frames by cluster
  list_per_cluster <- split(top_All_6cl, top_All_6cl$cluster)
  list_per_cluster2 <- split(top10, top10$cluster)
  # Extract the gene column from each data frame in the list
  genes_per_cluster_6cl <- lapply(list_per_cluster, function(df) df$gene)
  genes_per_cluster_2 <- lapply(list_per_cluster2, function(df) df$gene)
  
  
  # Function to calculate Jaccard Index
  jaccard_index <- function(set1, set2) {
    if (length(set1) == 0 || length(set2) == 0) {
      return(0)
    }
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    return(intersection / union)
  }
  
  # Assuming genes_per_cluster_6cl and genes_per_cluster_2 are your lists of lists
  # Calculate Jaccard Index for each combination of clusters between the two lists of lists
  jaccard_results <- matrix(0, nrow = length(genes_per_cluster_6cl), ncol = length(genes_per_cluster_2))
  rownames(jaccard_results) <- names(genes_per_cluster_6cl)
  colnames(jaccard_results) <- names(genes_per_cluster_2)
  
  for (i in seq_along(genes_per_cluster_6cl)) {
    for (j in seq_along(genes_per_cluster_2)) {
      jaccard_results[i, j] <- jaccard_index(genes_per_cluster_6cl[[i]], genes_per_cluster_2[[j]])
    }
  }
  
  # Convert the matrix to a long format suitable for ggplot
  jaccard_data_long <- as.data.frame(as.table(jaccard_results))
  
  # Renaming the columns for clarity
  colnames(jaccard_data_long) <- c("Cluster_6cl", "Cluster_2", "OverlapIndex")
  
  # Create the heatmap
  heat_jackard_6cl =    ggplot(jaccard_data_long, aes(x = Cluster_6cl, y = Cluster_2, fill = OverlapIndex)) +
    geom_tile(color="white") + # Use geom_tile() for heatmap squares
    scale_fill_viridis() + 
    theme_minimal() + # Use a minimal theme for cleaner appearance
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
    labs(fill = "jaccard\nIndex", x = NULL, y = NULL, title = "jaccard Index Heatmap") # Set labels and title
  
  
  ### OVERLAP ###
  # Function to calculate Overlap Index
  overlap_index <- function(set1, set2) {
    if (length(set1) == 0 || length(set2) == 0) {
      return(0)
    }
    intersection_size <- length(intersect(set1, set2))
    min_size <- min(length(set1), length(set2))
    return(intersection_size / min_size)
  }
  
  # Assuming genes_per_cluster_6cl and genes_per_cluster_2 are your lists of lists
  
  # Calculate Overlap Index for each combination of clusters between the two lists of lists
  overlap_results <- matrix(0, nrow = length(genes_per_cluster_6cl), ncol = length(genes_per_cluster_2))
  rownames(overlap_results) <- names(genes_per_cluster_6cl)
  colnames(overlap_results) <- names(genes_per_cluster_2)
  
  for (i in seq_along(genes_per_cluster_6cl)) {
    for (j in seq_along(genes_per_cluster_2)) {
      overlap_results[i, j] <- overlap_index(genes_per_cluster_6cl[[i]], genes_per_cluster_2[[j]])
    }
  }
  
  
  # Convert the matrix to a long format suitable for ggplot
  overlap_data_long <- as.data.frame(as.table(overlap_results))
  
  # Renaming the columns for clarity
  colnames(overlap_data_long) <- c("Cluster_6cl", "Cluster_2", "OverlapIndex")
  
  # Create the heatmap
  heat_overlap_6cl =    ggplot(overlap_data_long, aes(x = Cluster_6cl, y = Cluster_2, fill = OverlapIndex)) +
      geom_tile(color="white") + # Use geom_tile() for heatmap squares
      scale_fill_viridis() + 
    theme_minimal() + # Use a minimal theme for cleaner appearance
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
      labs(fill = "Overlap\nIndex", x = NULL, y = NULL, title = "Overlap Index Heatmap") # Set labels and title
    
  
  
  
  
  ######### Jackard and overlap Index for 3 clusters #############
  # Split the data frame into a list of data frames by cluster
  list_per_cluster <- split(top_All, top_All$cluster)
  list_per_cluster2 <- split(top10, top10$cluster)
  # Extract the gene column from each data frame in the list
  genes_per_cluster_3cl <- lapply(list_per_cluster, function(df) df$gene)
  genes_per_cluster_2 <- lapply(list_per_cluster2, function(df) df$gene)
  
  # Assuming genes_per_cluster_3cl and genes_per_cluster_2 are your lists of lists
  # Calculate Jaccard Index for each combination of clusters between the two lists of lists
  jaccard_results <- matrix(0, nrow = length(genes_per_cluster_3cl), ncol = length(genes_per_cluster_2))
  rownames(jaccard_results) <- names(genes_per_cluster_3cl)
  colnames(jaccard_results) <- names(genes_per_cluster_2)
  
  for (i in seq_along(genes_per_cluster_3cl)) {
    for (j in seq_along(genes_per_cluster_2)) {
      jaccard_results[i, j] <- jaccard_index(genes_per_cluster_3cl[[i]], genes_per_cluster_2[[j]])
    }
  }
  
  # Convert the matrix to a long format suitable for ggplot
  jaccard_data_long2 <- as.data.frame(as.table(jaccard_results))
  
  # Renaming the columns for clarity
  colnames(jaccard_data_long2) <- c("Cluster_3cl", "Cluster_2", "OverlapIndex")
  
  # Create the heatmap
  heat_jackard_3cl =    ggplot(jaccard_data_long2, aes(x = Cluster_3cl, y = Cluster_2, fill = OverlapIndex)) +
    geom_tile(color="white") + # Use geom_tile() for heatmap squares
    scale_fill_viridis() + 
    theme_minimal() + # Use a minimal theme for cleaner appearance
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
    labs(fill = "jaccard\nIndex", x = NULL, y = NULL, title = "jaccard Index Heatmap 3pops") # Set labels and title
  
  
  
  # Calculate Overlap Index for each combination of clusters between the two lists of lists
  overlap_results <- matrix(0, nrow = length(genes_per_cluster_3cl), ncol = length(genes_per_cluster_2))
  rownames(overlap_results) <- names(genes_per_cluster_3cl)
  colnames(overlap_results) <- names(genes_per_cluster_2)
  
  for (i in seq_along(genes_per_cluster_3cl)) {
    for (j in seq_along(genes_per_cluster_2)) {
      overlap_results[i, j] <- overlap_index(genes_per_cluster_3cl[[i]], genes_per_cluster_2[[j]])
    }
  }
  
  
  # Convert the matrix to a long format suitable for ggplot
  overlap_data_long2 <- as.data.frame(as.table(overlap_results))
  
  # Renaming the columns for clarity
  colnames(overlap_data_long2) <- c("Cluster_3cl", "Cluster_2", "OverlapIndex")
  
  # Create the heatmap
  heat_overlap_3cl =    ggplot(overlap_data_long2, aes(x = Cluster_3cl, y = Cluster_2, fill = OverlapIndex)) +
    geom_tile(color="white") + # Use geom_tile() for heatmap squares
    scale_fill_viridis() + 
    theme_minimal() + # Use a minimal theme for cleaner appearance
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
    labs(fill = "Overlap\nIndex", x = NULL, y = NULL, title = "Overlap Index Heatmap 3pops") # Set labels and title
  
  
  if (DO_PRINT_JACKARD==TRUE){
  cat(" \n \n")
  print(heat_jackard_3cl)
  cat(" \n \n")
  }
  
  if (DO_PRINT_OVERLAP==TRUE){
  cat(" \n \n")
  print(heat_overlap_3cl)
  cat(" \n \n")
  }
  
  
  
  palette_i = hue_pal()(length(levels(object@active.ident)))
  names(palette_i) = levels(object@active.ident)
  
  cat(" \n \n")
  cat("#### Scoring with 3 pops NK1 , NK2 and NK3")
  cat(" \n \n")
  
  
  
  p2= FeaturePlot(object, features="NK11" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NK1")   & NoAxes()
  p3= FeaturePlot(object, features="NK21" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NK2")   & NoAxes()
  p4=FeaturePlot(object, features="NK31" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  & ggtitle("NK3")   & NoAxes()
  
  
  cat(" \n \n")
  print((p1  + p2 + p3 + p4))
  cat(" \n \n")
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK1_2_3.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print((p1  + p2 + p3 + p4))
  dev.off()
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK1_2_3.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print((p1  + p2 + p3 + p4))
  dev.off()
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNKclust.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p1)
  dev.off()
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNKclust.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p1)
  dev.off()
  
  
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK1.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p2)
  dev.off()
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK1.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p2)
  dev.off()
  
  
  
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK2.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print( p3 )
  dev.off()
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK2.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p3)
  dev.off()
  
  
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK3.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p4)
  dev.off()
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNK3.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p4)
  dev.off()
  
  
  
  
  
  
  
  if (DO_LABEL_SAVE ==TRUE){
    
 
  
  p1_labeled = DimPlot(object, pt.size = PT_SIZE,  label=TRUE, repel = TRUE)  & ggtitle("Clust") & NoAxes()
    p1_big_groups= DimPlot(object, pt.size = PT_SIZE, group.by = "chief" ) & ggtitle("Clust") & NoAxes()
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNKclust_Labeled.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p1_labeled)
  dev.off()
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNKclust_Labeled.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p1_labeled)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNKclust_BigGroups.png") 
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p1_big_groups)
  dev.off()
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"FeaturePlotNKclust_BigGroups.pdf") 
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p1_big_groups)
  dev.off()
  
  
  
  }
  
  
  
  
  p5= VlnPlot(object, features="NK11" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i)  & ggtitle("NK1") & stat_summary(fun.data=data_summary,color="black")
  p6= VlnPlot(object, features="NK21" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i)  & ggtitle("NK2") & stat_summary(fun.data=data_summary,color="black")
  p7= VlnPlot(object, features="NK31" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i) & ggtitle("NK3") & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print((p5 ))
  cat(" \n \n")
  
  cat(" \n \n")
  print(( p6 ))
  cat(" \n \n")
  
  cat(" \n \n")
  print(( p7))
  cat(" \n \n")
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"VlnPlotNK1score.png")
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p5)
  dev.off()
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"VlnPlotNK1score.pdf")
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p5)
  dev.off()
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"VlnPlotNK2score.png")
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p6)
  dev.off()
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"VlnPlotNK2score.pdf")
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p6)
  dev.off()
  
  
  
  
  
  
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"VlnPlotNK3score.png")
  
  png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
  print(p7)
  dev.off()
  
  
  
  FILE_SAVE= paste0(PATH_ANALYSIS_OUTPUT_FIGUREPLOTS,"/", names(Co_List_Objects)[compteur_tour] ,"VlnPlotNK3score.pdf")
  
  pdf(file = FILE_SAVE, width = 15/2.54, height = 10/2.54)
  print(p7)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  p8= RidgePlot(object, features="NK11" , cols= palette_i, sort = SORT_RIDGE_PLOT )  & ggtitle("NK1")
  p9= RidgePlot(object, features="NK21" , cols= palette_i, sort = SORT_RIDGE_PLOT)  & ggtitle("NK2")
  p10= RidgePlot(object, features="NK31" , cols= palette_i, sort = SORT_RIDGE_PLOT) & ggtitle("NK3")
  
  cat(" \n \n")
  print((p8 ))
  cat(" \n \n")
  
  cat(" \n \n")
  print(( p9 ))
  cat(" \n \n")
  
  cat(" \n \n")
  print(( p10))
  cat(" \n \n")
  
  
  cat(" \n \n")
  cat("#### Scoring with 6 pops NK1A, NK1B, NK1C , NKint, NK2 and NK3")
  cat(" \n \n")
  
  p11 =DimPlot(object, pt.size = PT_SIZE) & ggtitle("Clust")
  p12= FeaturePlot(object, features="NK1A_6cl1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NK1A_6cl")   & NoAxes()
  p13= FeaturePlot(object, features="NK1B_6cl1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NK1B_6cl")   & NoAxes()
  p14= FeaturePlot(object, features="NK1C_6cl1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NK1C_6cl")   & NoAxes()
  p15= FeaturePlot(object, features="NKint_6cl1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NKint_6cl")   & NoAxes()
  p16= FeaturePlot(object, features="NK2_6cl1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & ggtitle("NK2_6cl")   & NoAxes()
  p17=FeaturePlot(object, features="NK3_6cl1" , pt.size = PT_SIZE) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  & ggtitle("NK3_6cl")   & NoAxes()
  
  
  if (DO_PRINT_JACKARD==TRUE){
    cat(" \n \n")
    print(heat_jackard_6cl)
    cat(" \n \n")
  }
  
  if (DO_PRINT_OVERLAP==TRUE){
    cat(" \n \n")
    print(heat_overlap_6cl)
    cat(" \n \n")
  }
  
  
  cat(" \n \n")
  print((p11  + p12 + p13 + p14))
  cat(" \n \n")
  
  cat(" \n \n")
  print((p11  + p15 + p16 + p17))
  cat(" \n \n")
  
  p25= VlnPlot(object, features="NK1A_6cl1" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT , cols= palette_i)  & ggtitle("NK1A_6cl") & stat_summary(fun.data=data_summary,color="black")
  p26= VlnPlot(object, features="NK1B_6cl1" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i)  & ggtitle("NK1B_6cl") & stat_summary(fun.data=data_summary,color="black")
  p27= VlnPlot(object, features="NK1C_6cl1" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i)  & ggtitle("NK1C_6cl") & stat_summary(fun.data=data_summary,color="black")
  p28= VlnPlot(object, features="NKint_6cl1" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i)  & ggtitle("NKint_6cl") & stat_summary(fun.data=data_summary,color="black")
  p29= VlnPlot(object, features="NK2_6cl1" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i)  & ggtitle("NK2_6cl") & stat_summary(fun.data=data_summary,color="black")
  p30= VlnPlot(object, features="NK3_6cl1" , pt.size = PT_SIZE_VLN, sort = SORT_VLN_PLOT, cols= palette_i) & ggtitle("NK3_6cl") & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p25 )
  cat(" \n \n")
  
  cat(" \n \n")
  print( p26 )
  cat(" \n \n")
  
  cat(" \n \n")
  print( p27)
  cat(" \n \n")
  
  cat(" \n \n")
  print( p28)
  cat(" \n \n")
  
  cat(" \n \n")
  print( p29)
  cat(" \n \n")
  
  cat(" \n \n")
  print( p30)
  cat(" \n \n")
  
  
  
  p35= RidgePlot(object, features="NK1A_6cl1"  , cols= palette_i, sort = SORT_RIDGE_PLOT)  & ggtitle("NK1A_6cl")
  p36= RidgePlot(object, features="NK1B_6cl1"  , cols= palette_i, sort = SORT_RIDGE_PLOT)  & ggtitle("NK1B_6cl")
  p37= RidgePlot(object, features="NK1C_6cl1"  , cols= palette_i, sort = SORT_RIDGE_PLOT)  & ggtitle("NK1C_6cl")
  p38= RidgePlot(object, features="NKint_6cl1" , cols= palette_i , sort = SORT_RIDGE_PLOT)  & ggtitle("NKint_6cl")
  p39= RidgePlot(object, features="NK2_6cl1"  , cols= palette_i, sort = SORT_RIDGE_PLOT)  & ggtitle("NK2_6cl")
  p40= RidgePlot(object, features="NK3_6cl1" , cols= palette_i, sort = SORT_RIDGE_PLOT )  & ggtitle("NK3_6cl")
  
  
  cat(" \n \n")
  print(p35 )
  cat(" \n \n")
  
  cat(" \n \n")
  print(p36 )
  cat(" \n \n")
  
  cat(" \n \n")
  print(p37)
  cat(" \n \n")
  
  cat(" \n \n")
  print(p38)
  cat(" \n \n")
  
  cat(" \n \n")
  print(p39)
  cat(" \n \n")
  
  cat(" \n \n")
  print(p40)
  cat(" \n \n")
  
  
  if(DO_MELSEN_Heatmap ==TRUE){
  df_For_Heatmap = object@meta.data[ ,c("NK11", "NK21", "NK31", "final_nomenclature")]
  
  aggregate_data <- aggregate(cbind(NK11, NK21, NK31) ~ final_nomenclature, data = df_For_Heatmap, mean)
  
  
  pheatmap(aggregate_data[,-1], 
           scale = "column", 
           annotation_names_row = TRUE, 
           annotation_names_col = TRUE, 
           annotation_row = data.frame(final_nomenclature = aggregate_data$final_nomenclature),
           show_rownames = TRUE, 
           show_colnames = TRUE,
           main = "Heatmap of NK11, NK21, and NK31 by Tissue Chief")
  }
  
}


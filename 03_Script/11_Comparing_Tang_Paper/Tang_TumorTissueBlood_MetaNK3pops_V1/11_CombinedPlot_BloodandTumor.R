#Try to make a combined plot for blood and Intra Tumor NK cells

#Load the query data
data.query_blood  = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_Blood_MetaNK3pops_V1/data.query_all_cells_BloodTang.rds")
data.query_tumor  = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_TumorTissueBlood_MetaNK3pops_V1/data.query_all_cells_Tumor_Normal.rds")
data.query_tumor = SplitObject(data.query_tumor , split.by = "meta_tissue")
data.query_tumor= data.query_tumor$Tumor


#Prepare Blood Data
data.query_blood$meta_histology= as.factor(data.query_blood$meta_histology)
data.query_blood$meta_histology_relabel  = data.query_blood$meta_histology
levels(data.query_blood$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(data.query_blood$meta_histology_relabel ))
data.query_blood$predicted.id =  as.factor(data.query_blood$predicted.id)
data.query_blood$SecondClust  = as.factor(data.query_blood$SecondClust)

#Prepare Tumor Data
data.query_tumor$meta_histology= as.factor(data.query_tumor$meta_histology)
data.query_tumor$meta_histology_relabel  = data.query_tumor$meta_histology
levels(data.query_tumor$meta_histology_relabel) = sub(".*\\((.*)\\).*", "\\1", levels(data.query_tumor$meta_histology_relabel ))
data.query_tumor$predicted.id =  as.factor(data.query_tumor$predicted.id)
data.query_tumor$SecondClust  = as.factor(data.query_tumor$SecondClust)



#Create proportion table for Tumor
Clusters_Proportions = prop.table( table( data.query_tumor$meta_histology , data.query_tumor$predicted.id) , margin = 1)
df = as.data.frame.matrix(Clusters_Proportions)
df$meta_histology <- rownames(df)
df_filtered <- df[!is.na(df$NK2), ]
df_ordered <- df_filtered[order(df_filtered$NK2), ]

  #Add a healthy donor row composed of NA
number_rows <- nrow(df_ordered)
df_ordered[number_rows+1,]= NA
row.names(df_ordered)[number_rows+1] = "Healthy donor"
df_ordered["Healthy donor","meta_histology"] = "Healthy donor"


ORDER_CANCER = rownames(df_ordered)
ORDER_CANCER_RELABEL = sub(".*\\((.*)\\).*", "\\1", ORDER_CANCER )



#Create proportion table for Blood
Clusters_Proportions2 = prop.table( table( data.query_blood$meta_histology , data.query_blood$predicted.id) , margin = 1)
df2 = as.data.frame.matrix(Clusters_Proportions2)
df2$meta_histology <- rownames(df2)


#Create a new df2_ordered for Blood with NA when not present
  #Initialize df2_ordered
df2_ordered <- data.frame(matrix(ncol = ncol(df2), nrow = nrow(df_ordered)))
colnames(df2_ordered) <- colnames(df2)


  # Set the meta_histology column in df2_ordered to match the order in df_ordered
df2_ordered$meta_histology <- df_ordered$meta_histology


  # Loop through df2_ordered to fill in values from df2 or NA
for (i in 1:nrow(df2_ordered)) {
  if (df2_ordered$meta_histology[i] %in% df2$meta_histology) {
    df2_ordered[i, ] <- df2[df2$meta_histology == df2_ordered$meta_histology[i], ]
  } else {
    df2_ordered[i,c("NK1","NK2","NK3") ] <- NA
  }
}

# Assign meta_histology as row names for easier readability
rownames(df2_ordered) <- df2_ordered$meta_histology


#OLD PLOT
p22 = ggplot(df2_ordered, aes(x=meta_histology, fill= predicted.id)) + geom_bar(position="fill")  + scale_fill_manual(values= palette3)+ ggtitle(unlist(paste0("Tang Only", title)) ) +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 



#Try based on CHATGPT suhggestion

# Convert the dataframes to long format
df_ordered_long <- df_ordered %>% 
  pivot_longer(cols = c(NK1, NK2, NK3), names_to = "Cluster", values_to = "Proportion") %>% 
  mutate(DataFrame = 'df_ordered')

df2_ordered_long <- df2_ordered %>% 
  pivot_longer(cols = c(NK1, NK2, NK3), names_to = "Cluster", values_to = "Proportion") %>% 
  mutate(DataFrame = 'df2_ordered')

# Combine the two dataframes
combined_df <- rbind(df_ordered_long, df2_ordered_long)

# Ensure meta_histology is a factor with levels in the same order as df_ordered
combined_df$meta_histology <- factor(combined_df$meta_histology, levels = df_ordered$meta_histology)

# Plot
p_combined = ggplot(combined_df, aes(x = meta_histology, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity",  na.rm = FALSE) +
  scale_fill_manual(values = c("NK1" = '#F8766D', "NK2" = '#8494FF', "NK3" = "#0CB702")) +
  facet_grid(DataFrame ~ ., scales = "free_y", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Meta Histology", y = "Proportion") +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 


FILE_SAVE= "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_Blood_MetaNK3pops_V1/Figures_Plots/CombinedPlot_TumorvsBlood/CombinedBarPlot.png"

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p_combined)
dev.off()



# Ensure meta_histology is a factor with levels in the same order as df_ordered
combined_df$meta_histology <- factor(combined_df$meta_histology, levels = df_ordered$meta_histology)
combined_df$meta_histology = sub(".*\\((.*)\\).*", "\\1", combined_df$meta_histology )
combined_df$meta_histology  = factor(combined_df$meta_histology, levels = sub(".*\\((.*)\\).*", "\\1", df_ordered$meta_histology ))

# Plot
p_combined = ggplot(combined_df, aes(x = meta_histology, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity",  na.rm = FALSE) +
  scale_fill_manual(values = c("NK1" = '#F8766D', "NK2" = '#8494FF', "NK3" = "#0CB702")) +
  facet_grid(DataFrame ~ ., scales = "free_y", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Meta Histology", y = "Proportion") +
  theme_classic()+
  theme(text = element_text(size = 7), axis.text.x = element_text(hjust=1, angle = 45) ) 


FILE_SAVE= "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/11_Comparing_Tang_Paper/Tang_Blood_MetaNK3pops_V1/Figures_Plots/CombinedPlot_TumorvsBlood/CombinedBarPlot_Abriviation.png"

png(file=FILE_SAVE, width = 15, height = 10,  units = "cm", res=600 )
print(p_combined)
dev.off()


#Fig7a
ggsave(
  "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/pdf_figures_NI/Fig7a.pdf",
  plot = p_combined,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 15,
  height = 10,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)




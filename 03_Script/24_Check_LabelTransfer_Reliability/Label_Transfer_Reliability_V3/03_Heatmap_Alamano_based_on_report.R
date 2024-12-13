
# 3. Melt the confusion matrix for ggplot

conf_matrix = matrix(c(8.6, 1.5, 86.4, 0.8, 86.8, 0, 90.7, 11.7, 13.6), nrow = 3, ncol = 3, byrow = TRUE)
rownames(conf_matrix) = c("NK3", "NK2", "NK1")
colnames(conf_matrix) = c("NK1", "NK2", "NK3")

# Melt the confusion matrix
conf_matrix_melted = melt(conf_matrix, varnames = c("True Label", "Predicted Label"))


plot_confu2 = ggplot(conf_matrix_melted, aes(x=`Predicted Label`, y=`True Label`, fill=value)) +
  geom_tile() +  # Create the heatmap tiles
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +  # Reverted "RdBu"
  geom_text(aes(label=sprintf("%.1f%%", value)), color="black", size=4) +  # Add text labels with percentage
  theme_minimal() +  # Use a minimal theme
  labs(x = "Predicted Label", y = "True Label", fill = "Percentage") +  # Label the axes and legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

print(plot_confu2)


png(file="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/24_Check_LabelTransfer_Reliability/Label_Transfer_Reliability_V3/HeatmapHD.png", width = 10, height = 7,  units = "cm", res=600 )
print(plot_confu2)
dev.off()


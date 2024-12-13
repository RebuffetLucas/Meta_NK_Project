#CHANTIER BOXPLOT
#INSERER EN LIGNE 510

#BoxPlot version
table(data_query_tissue$orig.ident, data_query_tissue$meta_histology)
Clusters_Proportions = prop.table( table( data_query_tissue$orig.ident , data_query_tissue$predicted.id) , margin = 1)
#Clusters_Proportions2 = t(Clusters_Proportions)
Clusters_Proportions = data.frame(Clusters_Proportions)

Clusters_Proportions2 = data_query_tissue@meta.data[,c("orig.ident", "meta_histology", "predicted.id")]
Clusters_Proportions2 = prop.table()


colnames(Clusters_Proportions) = c("orig.ident", "seurat_clusters", "Freq")

Clusters_Proportions$seurat_clusters = factor(Clusters_Proportions$seurat_clusters, levels= ORDER_CLUST_LEGEND)

p12 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_boxplot(outlier.shape =  NA, color = palette[ORDER_CLUST_LEGEND], lwd= 0.9, alpha= 0.1) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 3) +
  theme_classic() +  theme(text = element_text(size = 15))       



#Version BarPlot

p14 = Clusters_Proportions %>%
  ggplot(aes(x = seurat_clusters, y = Freq, group = seurat_clusters, color = orig.ident)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", fill= "white", color= palette[ORDER_CLUST_LEGEND], size=1.2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = "dodge", width = 0.25, size=0.6) +
  geom_jitter( position = position_jitter(0.4),  aes(fill = orig.ident) , colour= "black", shape = 21, inherit.aes = TRUE, size= 2) +
  theme_classic() +  theme(text = element_text(size = 15))       

p14




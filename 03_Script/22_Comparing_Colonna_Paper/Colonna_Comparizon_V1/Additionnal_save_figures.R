
p1 =DimPlot(object, pt.size = PT_SIZE, cols=palette_i, label = TRUE) & ggtitle("Clust")   & NoAxes()



pdf(file = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/22_Comparing_Colonna_paper/Colonna_Comparizon_V1/FiguresPlots/Fig6a2.pdf", width = 15/2.54, height = 10/2.54)
print(p1)
dev.off()



p2 =DimPlot(object, pt.size = PT_SIZE,  group.by = "chief") & ggtitle("Clust")  & NoAxes()



pdf(file = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/22_Comparing_Colonna_paper/Colonna_Comparizon_V1/FiguresPlots/Fig6a1.pdf", width = 15/2.54, height = 10/2.54)
print(p2)
dev.off()



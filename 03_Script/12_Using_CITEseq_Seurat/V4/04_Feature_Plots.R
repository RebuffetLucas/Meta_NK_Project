
FeaturePlot(Seurat_NK_Final, reduction= "wnn.umap", features= c("CD34", "CD45RA") , min.cutoff = 'q01', max.cutoff = 'q99') &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

FeaturePlot(Seurat_NK_Final, reduction= "wnn.umap", features= c("CD8", "CD2") , min.cutoff = 'q01', max.cutoff = 'q99') &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

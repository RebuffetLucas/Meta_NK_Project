#Trying to assign Cell cycle score


#Read Seurat Object
PBMC = readRDS(paste0(OBJECT_DIRECTORY_PATH, "/NK_All_15clusters.rds"))

#Score Cell Cycle
PBMC = CellCycleScoring(PBMC, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

#Vizualize
ggplot2::ggplot(PBMC@meta.data, aes(x=S.Score, y= G2M.Score)) + ggplot2::geom_point(aes(color=Phase)) + theme(axis.text.x = element_text(angle = 90))

#Re-assign
G1cells = Cells(subset(PBMC, subset= S.Score < 0.1 & G2M.Score < 0.2 ))
Scells = Cells(subset(PBMC, subset= S.Score > 0.1 & G2M.Score < 0.2 ))
G2Mcells = Cells(subset(PBMC, subset=  G2M.Score > 0.2 ))


#Modify naming in Seurat Object
PBMC$Phase[G1cells] = "G1"
PBMC$Phase[Scells] = "S"
PBMC$Phase[G2Mcells] = "G2M"

#Visalisation
ggplot2::ggplot(PBMC@meta.data, aes(x=S.Score, y= G2M.Score)) + ggplot2::geom_point(aes(color=Phase)) + theme(axis.text.x = element_text(angle = 90))

DimPlot(PBMC, group.by = "Phase" )



#Merging All the datasets together and running the different batch corrections and Obtaining the first clustering

#OpenAll

source("/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/Meta_NK_V5/03_Script/14_SCTonOurData/00A_Global_Dependencies.R")
source("/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/Meta_NK_V5/03_Script/14_SCTonOurData/00B_Analysis_Params.R")

#BiocManager::install("glmGamPoi")
#ForV2:
NoIPH = readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/00_RawData/Seurat_Obj_other_samples/Sent_By_Janine/merged_and_downsampled.rds")
NoIPH= subset(NoIPH, subset= Chemistry == "V2" )


#Merge them
Merged_Seurat= NoIPH
Merged_Seurat$orig.ident = droplevels(Merged_Seurat$orig.ident)

#SCTransform
sce_Merged.list <- SplitObject(Merged_Seurat, split.by = "orig.ident")
sce_Merged.list <- lapply(X = sce_Merged.list, FUN = SCTransform, method = "glmGamPoi")

#Integration
features <- SelectIntegrationFeatures(object.list = sce_Merged.list, nfeatures = 3000)
sce_Merged.list <- PrepSCTIntegration(object.list = sce_Merged.list, anchor.features = features)
sce_Merged.list <- lapply(X = sce_Merged.list, FUN = RunPCA, features = features)


immune.anchors <- FindIntegrationAnchors(object.list = sce_Merged.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)


immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)

immune.combined.sct <- FindNeighbors(immune.combined.sct, dims = 1:30, verbose = FALSE)

immune.combined.sct  =readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/16_SCTonOurDAta/SCT_rPCA_V2only.rds")

immune.combined.sct$seurat_clusters = NULL

immune.combined.sct <- FindClusters(immune.combined.sct, verbose = FALSE,  resolution= FINDCLUSTERS_RESOLUTION ,random.seed = SEED)


saveRDS(immune.combined.sct , file= "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/16_SCTonOurDAta/SCT_rPCA_V2only.rds")




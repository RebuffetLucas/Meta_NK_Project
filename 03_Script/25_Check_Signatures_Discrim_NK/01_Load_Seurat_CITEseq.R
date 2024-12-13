#This scripts makes sure that the 13 genes from Crinier 2018 discriminates well the NK populations and that the signature NK1, NK2 and NK3 is suited for NK only

#Load the data
library(SeuratDisk)

reference <- LoadH5Seurat("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/Seurat_CITEseq/pbmc_multimodal.h5seurat")
reference = SetIdent(reference,  value= "celltype.l1")

#Have a look at the clusters compo
table(reference$celltype.l1)
table(reference$celltype.l2)
table(reference$celltype.l3)

#Have a look at the metadata
table(reference$donor)
table(reference$time)
table(reference$lane)


#Use another object
PBMC = reference


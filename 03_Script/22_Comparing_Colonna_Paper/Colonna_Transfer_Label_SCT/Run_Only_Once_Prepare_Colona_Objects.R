#Preparing_Colonna's data

###Loading the data

#ILC Blood 
Co_ILCs_Blood_1a = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Blood_Fig1a.rds")
#DimPlot(Co_ILCs_Blood_1a, group.by = "final_nomenclature" )
#Co_ILCs_Blood_1a = SetIdent(Co_ILCs_Blood_1a, value = "final_nomenclature")

Co_ILCs_Blood_1a$meta_patientID = Co_ILCs_Blood_1a$orig.ident
Co_ILCs_Blood_1a$meta_histology = "Blood"
Co_ILCs_Blood_1a$meta_tissue = "Blood"
Co_ILCs_Blood_1a$celltype = Co_ILCs_Blood_1a$final_nomenclature
Co_ILCs_Blood_1a$Majortype = Co_ILCs_Blood_1a$final_nomenclature
Co_ILCs_Blood_1a$datasets = Co_ILCs_Blood_1a$orig.ident



#NK Blood
Co_NK_Blood_1f = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Blood_Fig1f.rds")
#DimPlot(Co_NK_Blood_1f, group.by = "final_nomenclature")
#Co_NK_Blood_1f = SetIdent(Co_NK_Blood_1f, value = "final_nomenclature")

Co_NK_Blood_1f$meta_patientID = "Co_single_NK_Blood"
Co_NK_Blood_1f$meta_histology = "Blood"
Co_NK_Blood_1f$meta_tissue = "Blood"
Co_NK_Blood_1f$celltype = Co_NK_Blood_1f$final_nomenclature
Co_NK_Blood_1f$Majortype = Co_NK_Blood_1f$final_nomenclature
Co_NK_Blood_1f$datasets = "Co_single_NK_Blood"



#NK Tonsil
Co_NK_Tonsil_2d = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Tonsil_Fig2d.rds")
#DimPlot(Co_NK_Tonsil_2d, group.by = "final_nomenclature")
#Co_NK_Tonsil_2d = SetIdent(Co_NK_Tonsil_2d, value = "final_nomenclature")

Co_NK_Tonsil_2d = SetIdent(Co_NK_Tonsil_2d, value = "final_nomenclature")
Co_NK_Tonsil_2d = subset(Co_NK_Tonsil_2d , idents= "8 -ILC3" , invert=TRUE )
Co_NK_Tonsil_2d = SCTransform(Co_NK_Tonsil_2d, verbose = FALSE)
Co_NK_Tonsil_2d$final_nomenclature = droplevels(Co_NK_Tonsil_2d$final_nomenclature)

Co_NK_Tonsil_2d$meta_patientID = as.factor(Co_NK_Tonsil_2d$orig.ident)
levels(Co_NK_Tonsil_2d$meta_patientID) = c("Tonsil_Other", "Tonsil_1"  ,    "Tonsil_2"     , "Tonsil_3"  ,    "Tonsil_4"  )
Co_NK_Tonsil_2d$meta_histology = "Tonsil"
Co_NK_Tonsil_2d$meta_tissue = "Tonsil"
Co_NK_Tonsil_2d$celltype = Co_NK_Tonsil_2d$final_nomenclature
Co_NK_Tonsil_2d$Majortype = Co_NK_Tonsil_2d$final_nomenclature
Co_NK_Tonsil_2d$datasets = Co_NK_Tonsil_2d$meta_patientID





#NK Lung
Co_NK_Lung_3a = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/Lung_Fig3a.rds")
#DimPlot(Co_NK_Lung_3a , group.by = "final_nomenclature")
#Co_NK_Lung_3a = SetIdent(Co_NK_Lung_3a, value = "final_nomenclature")
Co_NK_Lung_3a$final_nomenclature = as.factor(Co_NK_Lung_3a$final_nomenclature)
Co_NK_Lung_3a = SetIdent(Co_NK_Lung_3a, value = "final_nomenclature")
Co_NK_Lung_3a = subset(Co_NK_Lung_3a, idents = c("ILC1",  "ILC2_ILC3" ), invert = TRUE)
Co_NK_Lung_3a=SCTransform(Co_NK_Lung_3a,verbose = FALSE)
Co_NK_Lung_3a$final_nomenclature =droplevels(Co_NK_Lung_3a$final_nomenclature )

Co_NK_Lung_3a$meta_patientID = as.factor(Co_NK_Lung_3a$orig.ident)
Co_NK_Lung_3a$meta_histology = "Lung"
Co_NK_Lung_3a$meta_tissue = "Lung"
Co_NK_Lung_3a$celltype = Co_NK_Lung_3a$final_nomenclature
Co_NK_Lung_3a$Majortype = Co_NK_Lung_3a$final_nomenclature
Co_NK_Lung_3a$datasets = Co_NK_Lung_3a$meta_patientID




#NK Intestines
Co_ILCs_Intestines_4k = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/IEL_LPL_Fig4k.rds")
#DimPlot(Co_ILCs_Intestines_4k , group.by = "final_nomenclature")
Co_ILCs_Intestines_4k = SetIdent(Co_ILCs_Intestines_4k, value = "final_nomenclature")
levels(Co_ILCs_Intestines_4k$final_nomenclature)
Co_NK_Intestines_4k = subset(Co_ILCs_Intestines_4k , idents = c("ILC3"  ,  "ILC1" ,  "ILC3p",     "ILC3-ILC1" ) , invert=TRUE)
Co_NK_Intestines_4k= SCTransform(Co_NK_Intestines_4k,verbose = FALSE)
Co_NK_Intestines_4k$final_nomenclature = droplevels( Co_NK_Intestines_4k$final_nomenclature)
#DimPlot(Co_ILCs_Intestines_4k_NK_Only)

Co_NK_Intestines_4k$meta_patientID = as.factor(Co_NK_Intestines_4k$orig.ident)
Co_NK_Intestines_4k$meta_histology = "Intestines"
Co_NK_Intestines_4k$meta_tissue = "Intestines"
Co_NK_Intestines_4k$celltype = Co_NK_Intestines_4k$final_nomenclature
Co_NK_Intestines_4k$Majortype = Co_NK_Intestines_4k$final_nomenclature
Co_NK_Intestines_4k$datasets = Co_NK_Intestines_4k$meta_patientID



#ILCs all
Co_ILCs_All_Tissues_5a = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Colonna_Data/AllTissues_Fig5a_b_c.rds")
Co_ILCs_All_Tissues_5a$tissue_chief = as.factor(Co_ILCs_All_Tissues_5a$tissue_chief)
Co_ILCs_All_Tissues_5a = SetIdent(Co_ILCs_All_Tissues_5a, value = "tissue_chief")
Co_NK_All_Tissues_5a = subset( Co_ILCs_All_Tissues_5a, idents= c( "pbmc_ilcp" , "tonsil_ilcp"), invert= TRUE)
Co_NK_All_Tissues_5a= SCTransform(Co_NK_All_Tissues_5a,verbose = FALSE)
Co_NK_All_Tissues_5a$tissue_chief = droplevels(Co_NK_All_Tissues_5a$tissue_chief)


Co_NK_All_Tissues_5a$orig.ident = as.factor(Co_NK_All_Tissues_5a$orig.ident)
levels(Co_NK_All_Tissues_5a$orig.ident) = c("Lung_1_N"   ,   "Lung_2_F" ,     "Lung_3_N"  ,    "Lung_4_F"  ,    "N-1425-2_IEL",  "Tonsil_Other" ,"Tonsil_1"   ,   "Tonsil_2" ,    "Tonsil_3" ,     "Tonsil_4" )

Co_NK_All_Tissues_5a$meta_patientID = as.factor(Co_NK_All_Tissues_5a$orig.ident)
Co_NK_All_Tissues_5a$meta_histology = Co_NK_All_Tissues_5a$tissue_chief
Co_NK_All_Tissues_5a$meta_tissue = as.factor(Co_NK_All_Tissues_5a$tissue)
Co_NK_All_Tissues_5a$celltype = Co_NK_All_Tissues_5a$tissue_chief
Co_NK_All_Tissues_5a$Majortype = Co_NK_All_Tissues_5a$chief
Co_NK_All_Tissues_5a$datasets = Co_NK_All_Tissues_5a$orig.ident




#DimPlot(Co_ILCs_All_Tissues_5a, group.by = "chief")
#DimPlot(Co_ILCs_All_Tissues_5a, group.by = "tissue")
#DimPlot(Co_ILCs_All_Tissues_5a, group.by = "tissue_chief")


#Puting them all in a gene list
Co_List_Objects = list(Co_NK_Blood_1f , Co_NK_Tonsil_2d, Co_NK_Lung_3a, Co_NK_Intestines_4k)
names(Co_List_Objects) = c( "Co_NK_Blood_1f" , "Co_NK_Tonsil_2d", "Co_NK_Lung_3a", "Co_NK_Intestines_4k")

Co_List_Objects_Merged=  merge(Co_List_Objects[names(Co_List_Objects) ][[1]] , y = Co_List_Objects[names(Co_List_Objects) ][-1]) #Merge all the query together

DefaultAssay( Co_List_Objects_Merged ) = "RNA" 
saveRDS(Co_List_Objects_Merged, file = paste0(PATH_ANALYSIS_OUTPUT,"/Seurat_Objects/NKcellsfrom4tissues.rds"))



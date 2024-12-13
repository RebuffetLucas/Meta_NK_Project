###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis

######## CONSTANTS ON THE ANALYSIS STEP IDENTITY
######## Modify the ANALYSIS_STEP_NAME and LITERAL_TITLE variable

# This is the name of the analysis step. Use the same name as the folder
# name of the analysis step
# Example : ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity"
ANALYSIS_STEP_NAME = "24_Check_LabelTransfer_Reliability"
SUB_ANALYSIS_STEP_NAME= "Label_Transfer_Reliability_V1"


EXPERIMENT_NAME = SUB_ANALYSIS_STEP_NAME

# This is the literal title of the analysis step. It will be shown at the beginning
# of the HTML report
# Example : LITERAL_TITLE = "Quality Control, normalization and clustering"
LITERAL_TITLE = "Here we try to assess to what extent the label classic label transfer of Seurat is reliable we split MetaNK in 80% ref vs 20 % query and use the 4 pops of NK cells classification (NK1 , Nkint, NK2,  NK3) We train it by taking 500 cells per sample"

# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, SUB_ANALYSIS_STEP_NAME)


######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add here all the constants required for the analysis

#Analysis params
#Analysis params

SUBSET_DOWNSAMPLE_QUERY = TRUE #Whether downsample Tang Data
N_DOWNSAMPLE_QUERY_PER_IDENT = 500 #Number of cells to downasample randomly


#Building the reference
DATA_SPLIT =  "split" #split is paste0(Projet , orig.ident) avec projet = query / ref

#REFERENCE_FOR_INTEGRATION =  c("Dataset1", "Dataset2", "Dataset3", "Dataset4")

REFERENCE_FOR_INTEGRATION =  c("CMVneg1_donorA", "CMVneg2_donorB", "CMVpos2_donorC", "CMVpos3_donorD", "CMVpos4_donorE", "GSM3377678" ,"GSM3738542", "GSM3738543", "GSM5584154"  ,"GSM5584155", "GSM5584156_1", "GSM5584156_2", "GSM5584156_3")  

                               

#Defining variable features
#VARIABLE_FEATURES_METHOD = "FORCE_VF1" #FORCE_VF1 = take the intersection of VF1 and rownames as variable features
#"VARIABLE_FEATURE_BASIC" = Take classical 2000 VF within select integrations
VARIABLE_FEATURES_METHOD  = "VARIABLE_FEATURE_INTEGRATION_FEATURES" #=  Select Integration features 
VARIABLE_FEATURE_VAR_TO_SPLIT_TO_FIND_FEATURES= "orig.ident"


### Harmony batch Correction
VARIABLE_TANG_FOR_BATCHCOR = "meta_histology" #Batchcor = orig.ident for Meta_NK and (datasets, etc...) for Tang
Harmo_VARIABLE_BATCH_COR = c("Batch_Cor") 

HARMO_USE_REF = TRUE #True or false , WHether to use a reference during Harmony Integration
HARMO_REF_TO_USE  = c("CMVneg1_donorA", "CMVneg2_donorB" ,"CMVpos2_donorC", "CMVpos3_donorD", "CMVpos4_donorE" ,"GSM3377678" ,    "GSM3738542",
                      "GSM3738543" ,    "GSM5584154" ,    "GSM5584155"  ,   "GSM5584156_1" ,  "GSM5584156_2",   "GSM5584156_3"  )




NBCORES = 8


#Module scores
NUMBER_TOP_SCORING = 20


#Colors
palette<-c('NK1C'='#F8766D','NKint'='#8494FF',
           'NK1A'='#0CB702',
           'NK1B'='#00BFC4','NK2'='#ED68ED',
           'NK3'='#ABA300')

palette2<-c('NK1'='#F8766D','NKint'='#8494FF',
           'NK2'='#ED68ED',
           'NK3'='#ABA300')



# Seed for pseudo-random numbers
SEED = 42;


# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE
DATA_SCALE        = TRUE
DATA_VARS_REGRESS = NULL  # c("nCount_RNA") for UMIs (NULL to ignore)

#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report


# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.5;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden
FINDCLUSTERS_DIMS = 1:30

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)

FINDMARKERS_MINPCT    = 0.2      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10      # Number of marker genes to show in report and tables (NULL for all)

FDRCUTOFF = 0.05


###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "04_RNA_Velocity/Velocyto_Secondstep/CMVnegDonorA"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "04_Velocyto_CMVnegDonorA"


PATH_LOOM= file.path( PATH_PROJECT_OUTPUT,"Velocyto_FirstStep/CMVnegdonorA/possorted_genome_bam_LBQK3.loom")
PATH_SEURAT_RDS=file.path( "/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1")



FDRCUTOFF=0.05


SAMPLE  = "CMVneg1_donorA"

#Table of samples: CMVneg1_donorA CMVneg2_donorB CMVpos2_donorC CMVpos3_donorD CMVpos4_donorE     GSM3377678     GSM3738542     GSM3738543     GSM5584154     GSM5584155   GSM5584156_1   GSM5584156_2   GSM5584156_3 







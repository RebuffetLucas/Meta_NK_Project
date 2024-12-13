###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General


GLOBAL_DESCRIPTION = "scRNAseq_MetaNK_V5"

SCIENTIFIC_GROUP = "SUEVLAB"
SCIENTIFIC_PROJECT_NAME = "Meta_NK5"
EXPERIMENT_PROJECT_NAME = "Velocyto_CMVnegdonorA"



#### Additional custom variables (for a specific analysis step, tools, ...)

#SAMPLE_NAME = "CpG"
#SAMPLE_ID = ""




#### Input / Output

# RAW data folder (used to import raw data files into project folder)
PATH_RAWDATA = file.path( "/mnt/DOSI", 
                          SCIENTIFIC_GROUP,
                          "DATA/Genomics/BIOINFO_RAWDATA/ILCgroup")

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/mnt/DOSI", 
                          SCIENTIFIC_GROUP,
                          "BIOINFO/BIOINFO_PROJECT",
                          SCIENTIFIC_PROJECT_NAME)
                          #EXPERIMENT_PROJECT_NAME

PATH_PROJECT_RAWDATA      = file.path( PATH_PROJECT, "00_Rawdata")
PATH_PROJECT_EXTERNALDATA = file.path( PATH_PROJECT, "01_Reference")
PATH_PROJECT_CONTAINERS   = file.path( PATH_PROJECT, "02_Container")
PATH_PROJECT_OUTPUT       = file.path( PATH_PROJECT, "05_Output")





# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",
                            EXPERIMENT_PROJECT_NAME)
                            #, "_")
                            #startTimeFileName, "_",
                            #paramsHash, "_")

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



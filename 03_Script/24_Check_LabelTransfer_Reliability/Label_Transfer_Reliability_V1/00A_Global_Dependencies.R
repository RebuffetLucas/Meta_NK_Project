# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries
library(Matrix)
library( digest)        # digest (hashing)
library( DT)            # datatable
library( forcats)       # fct_inorder (factors)
library( fs)            # path_sanitize
library( future)        # plan (multicore)
library( ggplot2)
#library( ggpubr)
library( ggrepel)
#library( grid)
#library( gridExtra)
library( htmltools)     # browsable
library( htmlwidgets)   # JS (for datatable)
library( iheatmapr)     # iheatmap
library( knitr)
#library( magrittr)
library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
#library( plyr)
#library( dplyr)
library( rmarkdown)
library( scales)        # hue_pal

# Single-cell technology
library( Seurat)

# Functional Enrichment analysis
library( biomaRt)
library( clusterProfiler)
library( org.Mm.eg.db)
library( org.Hs.eg.db)
#library( rrvgo)

#Bertrand
library(gghalves)
library(ggdist)
library(patchwork)
library(ComplexHeatmap)

#Me
library(dplyr)
library(plyr)
library(readxl)
library(RColorBrewer)


#For batch correction with Seurat
library(BiocManager)
#library(multtest)
#library(metap)

#Harmony
library(harmony)
library(batchelor)
library(limma)
#library(unix)
library(DESeq2)
library(rlist)
library(reshape2)
library(rlist)

library(SeuratDisk)



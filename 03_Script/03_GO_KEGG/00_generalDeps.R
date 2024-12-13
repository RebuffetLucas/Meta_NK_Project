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
library( ggpubr)
library( ggrepel)
library( grid)
library( gridExtra)
library( htmltools)     # browsable
library( htmlwidgets)   # JS (for datatable)
library( iheatmapr)     # iheatmap
library( knitr)
library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
library( rmarkdown)
library( scales)        # hue_pal

# Single-cell technology
library( Seurat)

# Functional Enrichment analysis
library( biomaRt)
library( clusterProfiler)
library( org.Mm.eg.db)
library( org.Hs.eg.db)
library( rrvgo)

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

#Sabrina

library(webr)
library(ggforce)
library(data.table)
library(umap)
library(qgraph)
library(iheatmapr)
library(stringr)
library(future.apply)
library(qvalue)
library(factoextra)
library(knitr)
library(DT)
library(plotly)
library(ade4)
library(Rtsne)
library(tidyverse)
library(stringr)
library(scales)
library(dendsort)
library(dendextend)
library(VennDiagram)
library(memisc)
library(reshape2)
library(zoo)
library(viridis)
library(corrplot)
library(cowplot)
library(fs)
library(ReactomePA)
library(msigdbr)

#ANALYSIS PARAM

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 60;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000

#### Filtering / Normalization

# Filters for loading seurat object
LOAD_MIN_CELLS     = 3    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 200  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = 851
FILTER_UMI_MAX     = 21686

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = 603
FILTER_FEATURE_MAX = 4128

# Cells with percentage of mitochondrial genes above threshold will be excluded
FILTER_MITOPCT_MAX = 9.94

# Cells with percentage of ribosomal genes below threshold will be excluded
FILTER_RIBOPCT_MIN = 11.5

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize"
DATA_NORM_SCALEFACTOR = 10000

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE
DATA_SCALE        = FALSE
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

# PCA parameters
PCA_NPC              = 50  # Default number of dimensions to use for PCA (see Seurat::RunPCA())

PCA_PLOTS_NBDIMS     = 3   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15  # Number of 'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30  # Number of dimensions to use from PCA results
DIMREDUC_USE_UMAP_NBDIMS = 30  # Number of dimensions to use from UMAP results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)

FINDMARKERS_MINPCT    = 0.1      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 5       # Number of marker genes to show in report and tables (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05




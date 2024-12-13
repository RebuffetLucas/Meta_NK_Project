
# install.packages("reticulate")
library(reticulate)

use_condaenv("r-velocity", required = TRUE)
scv <- import("scvelo")
scv$logging$print_version()

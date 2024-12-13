# ####################################################################
# This script launch the compilation of report
# ####################################################################

library( knitr)
library( rmarkdown)
library( funr)


# Boolean to activate overwriting security (change at your own risk !)
preventOverwrite = TRUE;

### Define working folder (contains R/Rmd file for current sample, parent contains global project files)
if( exists( "snakemake")){
  WORKING_DIR = snakemake@scriptdir
}else{
  WORKING_DIR = dirname( sys.script())
}


### Load parameters
# Define an environment that will contain parameters
paramsEnv = new.env();

# Load file defining global parameters
globalParamsFilePath = file.path( WORKING_DIR, "../globalParams.R");
if(file.exists(globalParamsFilePath)) {
  source( globalParamsFilePath, local = paramsEnv);
} else {
  warning("The file 'globalParamsFilePath.R' containing global parameters is missing.");
}

# Load file defining sample parameters
sampleParamsFilePath = file.path( WORKING_DIR, "../sampleParams.R");
if(file.exists(sampleParamsFilePath)) {
  source( sampleParamsFilePath, local = paramsEnv);
} else {
  warning("The file 'sampleParamsFilePath.R' containing sample parameters is missing.");
}

# Load file defining analysis parameters
analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");
if(file.exists(analysisParamsFilePath)) {
  source( analysisParamsFilePath, local = paramsEnv);
} else {
  warning("The file 'analysisParamsFilePath.R' containing analysis-specific parameters is missing.");
}


### Prevent overwriting

# Create output file name
reportOutputFilename = paste0( paramsEnv[["SCIENTIFIC_PROJECT_NAME"]], "_",
                               paramsEnv[["EXPERIMENT_PROJECT_NAME"]], "_",
                               paramsEnv[["ANALYSIS_STEP_NAME"]], ".html");

alreadyExists = file.exists( file.path( paramsEnv[["PATH_ANALYSIS_OUTPUT"]], reportOutputFilename));
if( alreadyExists && preventOverwrite) stop( paste( c( "Report file already exists:", reportOutputFilename)));


### Prevent result mixing when parallel rendering
# Create a copy of original rmd file with unique name.
# Prevents 'render' to overwrite 'md' file created from name of 'Rmd' file,
# and scrambled results when parallel-rendering 'Rmd' files with identical name.
# Specifying 'tmp' folders for 'intermediates_dir' and 'knit_root_dir' prevent
# usage of WORKING_DIR in reports.
# Final output 'html' name provided to render anyway (not based on 'Rmd' name).

# Create a unique temporary file name
rmdCopyName = tempfile( pattern = "tempRmdCopy_", tmpdir = WORKING_DIR, fileext = ".Rmd");
# Copy 'Rmd' file to it
stopifnot( all( file.copy( from = file.path( WORKING_DIR, paste0( "Report_", paramsEnv[["ANALYSIS_STEP_NAME"]], ".Rmd")),
                           to   = rmdCopyName)));


### Render the report using previously built environment (use link to Rmd file)
rmarkdown::render( input = rmdCopyName,
                   output_dir = paramsEnv[["PATH_ANALYSIS_OUTPUT"]],
                   output_file  = I( reportOutputFilename),
                   envir = paramsEnv,
                   quiet = FALSE)

# Remove temporary 'Rmd' copy
if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));

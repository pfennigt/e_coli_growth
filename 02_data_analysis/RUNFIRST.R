# Adjust Paths to  the local investigation path &
# Install all needed packages

# Define parameters ---------------------------

# Define paths
# When structure is unchanged only INVESTIGATIONPATH has to be changed!
INVESTIGATIONPATH <- "C:/Users/tpeng/Desktop/e_coli_growth" #change to local investigation path
STUDYFOLDER <- "01_studies"
ANALYSISFOLDER <- "02_data_analysis"
FIGUREFOLDER <- "03_figures"

# Define required packages
required_packages <- c(
  "platexpress",
  "reshape2",
  "pracma",
  "ggplot2",
  "growthmodels",
  "growthrates",
  "pspline",
  "tidyr",
  "data.table",
  "plyr"
)

# Adjust data analysis scripts ---------------------------

# Find scripts
ANALYSISPATH <- file.path(INVESTIGATIONPATH, ANALYSISFOLDER)
analysis_files <- dir(ANALYSISPATH)

analysis_scripts <- grep("[0-9]{8}_analysis_SDY", analysis_files)
analysis_scripts <- analysis_files[analysis_scripts]

# Iterate adjustments through scripts
for (script_file in analysis_scripts) {
  
  # Read script to be adjusted
  current_script <- readLines(file.path(ANALYSISPATH, script_file))
  
  # Adjust paths ---------------------------
  
  # Get range of lines defining the paths
  path_range <- grep("# Define paths|# Set working directory", current_script)
  path_range[1] <- path_range[1]+1
  path_range[2] <- path_range[2]-2
  
  # Adjust paths withing these lines
  current_script[path_range[1]:path_range[2]] <- sprintf(
    "%s <- \"%s\"",
    c(
      "INVESTIGATIONPATH",
      "STUDYFOLDER",
      "ANALYSISFOLDER",
      "FIGUREFOLDER"
    ),
    c(
      INVESTIGATIONPATH,
      STUDYFOLDER,
      ANALYSISFOLDER,
      FIGUREFOLDER
    )
  )
  
  # Write adjusted script
  writeLines(
    text = current_script,
    con = file.path(ANALYSISPATH, script_file)
  )
}


# Install required packages ---------------------------

# Load devtools
require(devtools)

# Find the required packages not installed
installed_packages <- row.names(installed.packages())

new_packages <- required_packages[! required_packages %in% installed_packages]

# Install the non-installed required packages
for (package in new_packages) {

  # Thomas Petzoldt's `growthrates` package, development version
  if (package == "growthrates") {
    devtools::install_github('tpetzoldt/growthrates')
    
  # In-house package to parse and analyze platereader data
  } else if (package == "platexpress") {
    devtools::install_github('raim/platexpress')
  
  # In-house package adding growth models for use with `growthrates`
  } else if (package == "growthmodels") {
    devtools::install_git(
      'https://gitlab.com/raim/growthmodels.git',
      subdir = 'pkg',
      quiet = FALSE
    )
    
  # In-house package to segment growth curves into linear pieces
  } else if (package == "growthphases") {
    devtools::install_git(
      'https://gitlab.com/raim/growthphases.git',
      subdir = 'pkg',
      quiet = FALSE
    )
    
  } else {
    install.packages(package)
  }
}
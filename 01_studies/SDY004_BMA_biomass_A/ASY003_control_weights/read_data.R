# Read weight data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "drymass_weights"

# Define dataset study token
dat$study_id <- "SDY004"
dat$study_token <- "BMA"
dat$study_sub <- NA

# Define paths and locations
DATPATH <- "01_studies/SDY004_BMA_biomass_A/ASY003_control_weights"
data.file <- file.path(DATPATH, "20200605_control_weights_GA.csv")
meta.file <- file.path(DATPATH, "20200605_control_weights_GA_METADATA.csv")

# Read weight data
dat$data <- read.csv(data.file)

# Read metadata
dat$metadata = read.csv(meta.file)

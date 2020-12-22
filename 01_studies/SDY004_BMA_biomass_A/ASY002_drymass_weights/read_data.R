# Read weight data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define dataset study token
dat$study_id <- "SDY004"
dat$study_token <- "BMA"
dat$study_sub <- NA

# Define paths and locations
DATPATH <- "01_studies/SDY004_BMA_biomass_A/ASY002_drymass_weights"
data.file <- file.path(DATPATH, "20200605_drymass_weights_A.csv")
meta.file <- file.path(DATPATH, "20200605_drymass_weights_A_METADATA.csv")

# Read weight data
dat$data <- read.csv(data.file)

# Read metadata
dat$metadata = read.csv(meta.file)

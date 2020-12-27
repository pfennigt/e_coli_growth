# Read weight data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "drymass_weights"

# Define dataset study token
dat$study_id <- "SDY002"
dat$study_token <- "BMG"
dat$study_sub <- NA

# Define paths and locations
DATPATH <- "01_studies/SDY002_BMG_biomass_G/ASY002_drymass_weights"
data.file <- file.path(DATPATH, "20200303_drymass_weights_G.csv")
meta.file <- file.path(DATPATH, "20200303_drymass_weights_G_METADATA.csv")

# Read weight data
dat$data <- read.csv(data.file)

# Read metadata
dat$metadata = read.csv(meta.file)

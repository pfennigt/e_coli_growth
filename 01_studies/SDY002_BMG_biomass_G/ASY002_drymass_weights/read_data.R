# Read weight data
# Have working directory set to investigation folder

# Define paths and locations
DATPATH <- "01_studies/SDY002_BMG_biomass_G/ASY002_drymass_weights"
data.file <- file.path(DATPATH, "20200303_drymass_weights_G.csv")
data.file2 <- file.path(DATPATH, "20200303_drymass_weights_G_METADATA.csv")

# Read weight data
dat <- list()
dat$data <- read.csv(data.file)

dat$metadata = read.csv(data.file2)

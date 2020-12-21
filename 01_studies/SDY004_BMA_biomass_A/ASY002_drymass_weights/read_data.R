# Read weight data
# Have working directory set to investigation folder

# Define paths and locations
DATPATH <- "01_studies/SDY004_BMA_biomass_A/ASY002_drymass_weights"
data.file <- file.path(DATPATH, "20200605_drymass_weights_A.csv")
data.file2 <- file.path(DATPATH, "20200605_drymass_weights_A_METADATA.csv")

# Read weight data
dat <- list()
dat$data <- read.csv(data.file)

dat$metadata = read.csv(data.file2)

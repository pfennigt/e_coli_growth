# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "metabolite_measurements"

# Define dataset study token
dat$study_id <- "SDY010"
dat$study_token <- "MetabG"
dat$study_sub <- NA

# Define glucose as carbon source for plot descriptions
dat$CSource = "glucose"
dat$CSourceRead = "glc"

# Define paths and locations
DATPATH <- "01_studies/SDY010_MetabG_metabolites_G/ASY003_metabolite_glucose_glycogen_test"
data.file <- file.path(DATPATH, "20201016_metabolite_glucose_glycogen_test_best_gain.CSV")
meta.file <- file.path(DATPATH, "20201016_metabolite_glucose_glycogen_test_METADATA.csv")

# Read weight data
dat$data <- read.csv(data.file, skip = 14, header=T)

# Read metadata
dat$metadata <- read.csv(meta.file)
# Read BioLector data for glucose wells
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define dataset study token
dat$study_id <- "SDY002"
dat$study_token <- "BMG"
dat$study_sub <- NA

# Define glucose as carbon source for plot descriptions
dat$CSource = "glucose"
dat$CSourceRead = "glc"

# Define paths and locations
expid <- "27_Ecoli_2020_REDUCTION-1"
DATPATH <- "01_studies/SDY002_BMG_biomass_G/ASY001_BioLector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))

# Define used colors
Amounts = "amount"
PltColors = "color"
dat$bgCols = "royalblue3"

# Define analyzed time range
dat$xrng = c(0,12.5)

# Define parameters for the scatter peak search
dat$ignoreWells = c(paste0(rep(LETTERS[1:6], each=7), 1:7), "A8")
dat$peakTimeRange = NA

# Read biolector data
dat$data <- readExperiment(data.file,
                      type = "BioLectorPro",
                      time.conversion = 1/3600,
                      layout = layout.file, 
                      fields = c("strain","glc","ace","aTc"), 
                      afields = c("glc","ace","aTc"),
                      blank.id = "blank",
                      blank.data = c("Biomass","Riboflavine", "NADH - NADPH"),
                      group1 = dat$CSourceRead,
                      group2 = c(dat$CSourceRead,Amounts),
                      group2.color = PltColors
)

# Read plot annotation data
dat$anno <- NA
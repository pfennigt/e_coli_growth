# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define glucose as carbon source for plot descriptions
dat$CSource = "acetate"
dat$CSourceRead = "ace"

# Define paths and locations
expid <- "31_Ecoli_2020-1831_REDUCTION-1"
DATPATH <- "01_studies/SDY004_BMA_biomass_A/ASY001_biolector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
annotation.file <- file.path(DATPATH,"plot_annotations.csv")

# Define used colors
Amounts = "ace.amount"
PltColors = "ace.color"
dat$bgCol = "firebrick2"

# Define analyzed time range
dat$xrng = c(0,40)

# Define parameters for the scatter peak search
dat$ignoreWells = paste0(rep(LETTERS[1:6], each=7), 1:7)
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
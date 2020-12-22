# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define dataset study token
dat$study_id <- "SDY009"
dat$study_token <- "Z1WT"
dat$study_sub <- "WT"

# Define glucose as carbon source for plot descriptions
dat$CSource = "glucose"
dat$CSourceRead = "glc"

# Define paths and locations
expid <- "36_Ecoli_2020-1831_REDUCTION-1"
DATPATH <- "01_studies/SDY007_HighA_gradient_A/ASY001_biolector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))

# Define skipped wells (no WT data)
skipWellList = c(paste0(toupper(letters[1:3]), rep(1:8,3)))

# Define used colors
Amounts = "amount"
PltColors = "color"
dat$bgCols = colorRampPalette(c("white", "royalblue3"))(7)

# Define analyzed time range
dat$xrng = c(0,15)

# Define parameters for the scatter peak search
dat$ignoreWells = c(paste0(LETTERS[4:6], 1), "D8","E8","F8")
dat$peakTimeRange = c(3,10)

# Read biolector data
dat$data <- readExperiment(data.file,
                           type = "BioLectorPro",
                           time.conversion = 1/3600,
                           layout = layout.file, 
                           fields = c("strain","glc","ace","aTc"), 
                           afields = c("glc","ace","aTc"),
                           blank.id = "blank",
                           blank.data = c("Biomass","Riboflavine", "NADH - NADPH"),
                           skip.wells  =  skipWellList,
                           group1 = dat$CSourceRead,
                           group2 = c(dat$CSourceRead,Amounts),
                           group2.color = PltColors
)

# Read plot annotation data
dat$anno <- NA
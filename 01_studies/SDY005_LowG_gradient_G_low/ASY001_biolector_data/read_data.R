# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define dataset study token
dat$study_id <- "SDY005"
dat$study_token <- "LowG"
dat$study_sub <- NA

# Define glucose as carbon source for plot descriptions
dat$CSource = "glucose"
dat$CSourceRead = "glc"

# Define paths and locations
expid <- "32_Ecoli_2020-1831_REDUCTION-1"
DATPATH <- "01_studies/SDY005_LowG_gradient_G_low/ASY001_biolector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
annotation.file <- file.path(DATPATH,"plot_annotations.csv")

# Define used colors
Amounts = "amount"
PltColors = "color"
dat$bgCols = colorRampPalette(c("white", "royalblue3"))(15)

# Define analyzed time range
dat$xrng = c(0,15)

# Define parameters for the scatter peak search
dat$ignoreWells = c(paste0(rep(LETTERS[1:3], each=4), 1:4), "D8","E8","F8")
dat$peakTimeRange = c(5.5, 14.5)

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
dat$anno <- read.csv(annotation.file, header = T)
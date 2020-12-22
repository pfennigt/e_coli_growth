# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define dataset study token
dat$study_id <- "SDY006"
dat$study_token <- "InocG"
dat$study_sub <- NA

# Define glucose as carbon source for plot descriptions
dat$CSource = "glucose"
dat$CSourceRead = "glc"

# Define paths and locations
expid <- "33_Ecoli_2020-1831_REDUCTION-1"
DATPATH <- "01_studies/SDY006_InocG_inoc_gradient_G/ASY001_biolector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
annotation.file <- file.path(DATPATH,"plot_annotations.csv")

# Define used colors
Amounts = "inoc.amount"
PltColors = "inoc.color"
dat$bgCols = colorRampPalette(c("royalblue3", "green4"))(500)[c(20, 25, 32, 40, 50, 63, 79, 100, 126, 158, 199, 251, 316, 397, 500)]

# Define analyzed time range
dat$xrng = c(0,18.5)

# Define parameters for the scatter peak search
dat$ignoreWells = c("D8","E8","F8")
dat$peakTimeRange = c(5, 17)

# Read biolector data
dat$data <- readExperiment(data.file,
                           type = "BioLectorPro",
                           time.conversion = 1/3600,
                           layout = layout.file, 
                           fields = c("strain","glc","ace","aTc","inoc"), 
                           afields = c("glc","ace","aTc","inoc"),
                           blank.id = "blank",
                           blank.data = c("Biomass","Riboflavine", "NADH - NADPH"),
                           group1 = dat$CSourceRead,
                           group2 = c(dat$CSourceRead,Amounts),
                           group2.color = PltColors
)

# Read plot annotation data
dat$anno <- read.csv(annotation.file, header = T)
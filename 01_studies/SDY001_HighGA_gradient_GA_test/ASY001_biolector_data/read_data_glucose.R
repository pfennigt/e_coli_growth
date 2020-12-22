# Read BioLector data for glucose wells
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

# Define dataset study token
dat$study_id <- "SDY001"
dat$study_token <- "HighGA_glucose"
dat$study_sub <- "glucose"

# Define glucose as carbon source for plot descriptions
dat$CSource = "glucose"
dat$CSourceRead = "glc"

# Define paths and locations
expid <- "25_Ecoli_2020_REDUCTION-1"
DATPATH <- "01_studies/SDY001_HighGA_gradient_GA_test/ASY001_BioLector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
annotation.file <- file.path(DATPATH,"plot_annotations_glucose.csv")

# Define skipped wells (no acetate data)
skipWellList = c(paste0(toupper(letters[4:6]), rep(1:7,3)))

# Define used colors
Amounts = "amount"
PltColors = "color"
dat$bgCols = colorRampPalette(c("white", "royalblue3"))(7)

# Define analyzed time range
dat$xrng=c(0,12.5)

# Define parameters for the scatter peak search
dat$ignoreWells = paste0(rep(LETTERS[1:6], each=2), c(1,8))
dat$peakTimeRange = c(5.5, 12.5)

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
dat$anno <- read.csv(annotation.file, header = T)
# Read BioLector data for acetate wells
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "biolector_data"

# Define dataset study token
dat$study_id <- "SDY001"
dat$study_token <- "HighGA"
dat$study_sub <- "acetate"

# Define paths and locations
expid <- "25_Ecoli_2020_REDUCTION-1"
DATPATH <- "01_studies/SDY001_HighGA_gradient_GA_test/ASY001_BioLector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
annotation.file <- file.path(DATPATH,"plot_annotations_acetate.csv")

# Define experimental parameters (formatted to be parsed)
dat$carbon_source = "acetate"
dat$varied_factor = "acetate~concentration"
dat$varied_factor_unit = "C*'-'*mmol~l^-1"

# Define used BioLector parameters
dat$biolector_parameters = list(
  # specify the renaming of the raw data names
  measurements = c( 
    scatter="Biomass",
    ribof="Riboflavine",
    O2="DO(Pst3)",
    pH="pH(HP8)",
    NADH="NADH - NADPH"
  ),
  # units of the raw data measurements (formatted to be parsed)
  units = c( 
    "scatter" = "AU",
    "ribof" = "AU",
    "O2" = "'%'",
    "pH" = "' '",
    "NADH" = "AU"
  ),
  # descriptive text for the raw data measurements (formatted to be parsed)
  descriptors = c( 
    "scatter" = "backscatter",
    "ribof" = "riboflavin",
    "O2" = "O[2]~saturation",
    "pH" = "pH",
    "NADH" = "NADH"
  ),
  kLa = 230, # 1/h
  well_volume = 1 # ml
)

# Define skipped wells (no glucose data)
skipWellList = c(paste0(toupper(letters[1:3]), rep(1:7,3)))

# Define internal grouping
dat$layout_group = "ace"
dat$layout_amount = "ace.amount"
dat$layout_color = "ace.color"

# Define used colors
dat$bgCols = colorRampPalette(c("white", "firebrick2"))(7)

# Define analyzed time range
dat$xrng=c(0,37.5)

# Define parameters for the scatter peak search
dat$ignoreWells = paste0(rep(LETTERS[1:6], each=2), c(1,8))
dat$peakTimeRange = c(15, 25)

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
                      group1 = dat$layout_group,
                      group2 = c(dat$layout_group,dat$layout_amount),
                      group2.color = dat$layout_color
)

# Read plot annotation data
dat$anno <- read.csv(annotation.file, header = T)
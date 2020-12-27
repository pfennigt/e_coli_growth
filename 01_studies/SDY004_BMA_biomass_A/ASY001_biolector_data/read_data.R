# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "biolector_data"

# Define dataset study token
dat$study_id <- "SDY004"
dat$study_token <- "BMA"
dat$study_sub <- NA

# Define paths and locations
expid <- "31_Ecoli_2020-1831_REDUCTION-1"
DATPATH <- "01_studies/SDY004_BMA_biomass_A/ASY001_biolector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))

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
    "pH" = "",
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

# Define internal grouping
dat$layout_group = "ace"
dat$layout_amount = "ace.amount"
dat$layout_color = "ace.color"

# Define used colors
dat$bgCols = "firebrick2"

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
                           group1 = dat$layout_group,
                           group2 = c(dat$layout_group,dat$layout_amount),
                           group2.color = dat$layout_color
)

# Read plot annotation data
dat$anno <- NA
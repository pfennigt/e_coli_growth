# Read BioLector data for glucose wells
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "biolector_data"

# Define dataset study token
dat$study_id <- "SDY002"
dat$study_token <- "BMG"
dat$study_sub <- NA

# Define paths and locations
expid <- "27_Ecoli_2020_REDUCTION-1"
DATPATH <- "01_studies/SDY002_BMG_biomass_G/ASY001_BioLector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))

# Define experimental parameters (formatted to be parsed)
dat$carbon_source = "glucose"
dat$varied_factor = NA
dat$varied_factor_unit = NA

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
dat$layout_group = "glc"
dat$layout_amount = "amount"
dat$layout_color = "color"

# Define used colors
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
                      group1 = dat$layout_group,
                      group2 = c(dat$layout_group,dat$layout_amount),
                      group2.color = dat$layout_color
)

# Read plot annotation data
dat$anno <- NA
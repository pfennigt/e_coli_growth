# Read BioLector data
# Have working directory set to investigation folder

# Save data in a list
dat <- list()

#Define data type
dat$data_type = "biolector_data"

# Define dataset study token
dat$study_id <- "SDY006"
dat$study_token <- "InocG"
dat$study_sub <- NA

# Define paths and locations
expid <- "33_Ecoli_2020-1831_REDUCTION-1"
DATPATH <- "01_studies/SDY006_InocG_inoc_gradient_G/ASY001_biolector_data"
data.file <- file.path(DATPATH, paste0(expid,".csv"))
layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
annotation.file <- file.path(DATPATH,"plot_annotations.csv")

# Define experimental parameters (formatted to be parsed)
dat$carbon_source = "glucose"
dat$varied_factor = "inoculation~OD[600]"
dat$varied_factor_unit = "' '"

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

# Define internal grouping
dat$layout_group = "inoc"
dat$layout_amount = "inoc.amount"
dat$layout_color = "inoc.color"

# Define used colors
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
                           group1 = dat$layout_group,
                           group2 = c(dat$layout_group,dat$layout_amount),
                           group2.color = dat$layout_color
)

# Read plot annotation data
dat$anno <- read.csv(annotation.file, header = T)
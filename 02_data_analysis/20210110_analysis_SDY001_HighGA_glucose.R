# Initialize ---------------------------

# Load libraries
library(platexpress)
library(reshape2)
library(pracma)
library(ggplot2)
library(mugro)
library(growthrates)
library(pspline)
library(tidyr)
library(data.table)
library(plyr)

# Define analyzed study
study_id <- "SDY001"
study_sub <- "glucose"

# Define used drymass estimation
drymass_model <- "SDY002_BMG_drymass_estimation_twophase_ribof"

# Define paths
INVESTIGATIONPATH <- "C:/Users/tpeng/Desktop/e_coli_growth"
STUDYFOLDER <- "01_studies"
ANALYSISFOLDER <- "02_data_analysis"
FIGUREFOLDER <- "03_figures"

# Set working directory
setwd(INVESTIGATIONPATH)

# Get current date
current_date <- format(Sys.time(), "%Y%m%d")

# Initialize lists to store data and plots
dats <- list()
plts <- list()
plts$objects <- list()

# Load custom functions
custom_functions_file <- grep(
  pattern = "custom_functions",
  x = dir(ANALYSISFOLDER),
  value = T
)
source(file.path(ANALYSISFOLDER, custom_functions_file))


# Get study info ---------------------------

# Retrieve scurrent study from index
studies = read.csv(file.path(STUDYFOLDER, "studies_index.csv"))

current_study <- 
  studies$study_id == study_id & 
  if (is.na(study_sub)) {
    TRUE
  } else {
    studies$study_sub == study_sub
  }

current_study = studies[current_study, ]

# Set data path for current study
STUDYPATH <- file.path(STUDYFOLDER, current_study$study_path)


# Analysis of BioLector data  ---------------------------

# Define assay path
assay_id <- "ASY001"
assay_term <- sprintf("^%s[^\\.]*$", assay_id)

ASSAYFOLDER <- dir(STUDYPATH)[
  grep(assay_term, dir(STUDYPATH))
]
ASSAYPATH <- file.path(STUDYPATH,ASSAYFOLDER)

# Define read_data file
READFILE <- if (is.na(study_sub)) {
  "read_data.R"
} else {
  sprintf("read_data_%s.R", study_sub)
}

# Read data
source(file.path(ASSAYPATH, READFILE))

# Rename measurements
dat$data <- prettyData(dat$data, yids=dat$biolector_parameters$measurements)

# Adjust the layout
dat$data$layout[[dat$layout_amount]][dat$data$layout$strain=="blank"] <- NA
dat$data$layout$color[dat$data$layout$strain=="blank"] <- "#999999"

# Cut BioLector time range to predefined interesting time range
dat$data <- cutData(dat$data, xrng=dat$xrng)

# Get BioLector time points
times = dat$data$Time
nTime = length(times)

# Get the BioLector wells and their correct order
well_order = order(dat$data$layout$well)
wells = dat$data$layout$well[well_order]

# Get the amounts of the varied factor
factor_amount = dat$data$layout[[dat$layout_amount]][well_order]
factor_amount_uniq = sort(unique(factor_amount))

# Save in a list
dat$dataParameters <- list(
  times = times,
  nTime = nTime,
  well_order = well_order,
  wells = wells,
  factor_amount = factor_amount,
  factor_amount_uniq = factor_amount_uniq
)

# Define conversion factors
cpg <- 0.474 # g/g, gram carbon per gram biomass, Folsom&Carlson 2015
mpc <- 1/12 # mol/g, mol carbon per gram carbon
gps <- 0.603 # mg/ml/scatter, biomass per scatter

# Select non sampled wells for analysis
reference_wells <- 
  ! dat$data$layout$well %in% dat$ignoreWells & 
  dat$data$layout$strain != "blank"
reference_wells = wells[reference_wells[well_order]]

# Find drops in the scatter measurements (only time points after 2 h are considered)
scatter_drops <- findDrops(
  mat = getData(dat$data, "scatter")[, reference_wells],
  times = times,
  groups = factor_amount[wells %in% reference_wells],
  derivQuantThresh = 0.05,
  derivQuantThreshExpand = 0.1,
  inclTimeRange = 2,
  exclLast = T,
  doPlt = F
)

scatter_drops <- combineDrops(scatter_drops)

# Save peaks as a ggplot annotation
plts$objects[["scatter_drop_annotation"]] <- geom_rect(
  data = scatter_drops,
  aes(xmin = t1, xmax = t2),
  fill= "grey",
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.5
)

# Get values of ribof and scatter at the drop limits
scatter_drop_values = getDropValues(
  dat$data,
  IDs = c("scatter", "ribof"),
  combined_drops = scatter_drops,
  wells = reference_wells
)
names(scatter_drop_values)[3] <- "value_type"

# Save BioLector data for later use
dats$biolector <- dat


# Load drymass estimation  ---------------------------

drymass_lms_filename <- sprintf("%s.rds", drymass_model)
drymass_lms <- readRDS(
  file.path(
    ANALYSISFOLDER,
    "data_analysis_outputs",
    drymass_lms_filename
  )
)


# Continue analysis of BioLector data  ---------------------------

# Get the time of the first scatter drop in each well for dry mass estimation
scatter_drops_wells <- sapply(factor_amount, function(amount){
  if (is.na(amount) | amount == 0) {
    return(Inf)
    
  } else {
    res = scatter_drops[scatter_drops$groups == amount, "t1"][1]
    
    if (is.na(res)) {
      return(Inf)
      
    } else {
      return(res)
    }
  }
})


# Calculate the biomass concentration [gDW/l] and save
drymass_estimator_data <- getData(dats$biolector$data, drymass_lms$estimator)

estimated_drymass_g_l <- lapply(
  1:ncol(drymass_estimator_data),
  function(x){
    estimate_two_phases(
      values = drymass_estimator_data[,x],
      times = times,
      switchTime = scatter_drops_wells[x],
      lm1 = drymass_lms$lm1,
      lm2 = drymass_lms$lm2
    )
  }
)

estimated_drymass_g_l <- do.call(cbind, estimated_drymass_g_l)
dimnames(estimated_drymass_g_l) <- dimnames(drymass_estimator_data)

dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "drymass_g",
  dat = estimated_drymass_g_l
)

# Calculate the biomass amount [gDW]
estimated_drymass_g <- estimated_drymass_g_l * 0.001 #gDW

# Convert to C-mmol/l and save
estimated_drymass <- estimated_drymass_g_l * cpg * mpc * 1000 # C-mmol/l

dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "drymass",
  dat = estimated_drymass
)

# Adjust O2 saturation to max 100 %
O2 <- getData(dats$biolector$data, "O2")
O2 <- O2/max(O2)*100 # %

dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "O2",
  dat = O2,
  replace = T
)

# Calculate O2 concentration
O2_conc <- apply(O2, 2, O2PercToConc) # mmol/l

dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "O2_conc",
  dat = O2_conc
)

# Calculate H+ concentration [mol/l]
pH <- getData(dats$biolector$data,"pH")

HPlus <- 10^(-pH) # mol/l

dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "HPlus",
  dat = HPlus
)

# Add the new units and descriptors
dats$biolector$biolector_parameters$units <- c(
  dats$biolector$biolector_parameters$units,
  drymass_g = "gDW~l^-1",
  drymass = "C*'-'*mmol~l^-1",
  O2_conc = "mmol~l^-1",
  HPlus = "mol~l^-1",
  O2_consumption = "mmol~h^-1~gDW^-1"
)

dats$biolector$biolector_parameters$descriptors <- c(
  dats$biolector$biolector_parameters$descriptors,
  drymass_g = "dry~mass~conc.",
  drymass = "dry~mass~conc.",
  O2_conc = "O[2]~conc.",
  HPlus = "H^+~conc.",
  O2_consumption = "O[2]~consumption~rate"
)

# Calculate normalized measurements and normalized rates
measures <- names(dats$biolector$biolector_parameters$units)
measures <- measures[! measures %in% c("drymass", "drymass_g", "O2_consumption")]

for(meas in measures){
  
  # Get the data to be normalized
  measure_data <- getData(dats$biolector$data, meas)
  
  # Normalize with drymass in grams and save
  measure_norm <- measure_data / estimated_drymass_g #X/gDW
  
  dats$biolector$data <- addData(
    data = dats$biolector$data,
    ID = paste0(meas, "_norm"),
    dat = measure_norm
  )
  
  # Smooth the data with sm.spline
  measure_smooth <- smoothDat(
    varDat = measure_data,
    time = times,
    nOrder = 3
  )$smthObj
  
  # Predict the first time-derivative
  measure_deriv <- predictSMMat(measure_smooth, times, 1)
  
  # Normalize the time-derivative with drymass in grams and save
  measure_deriv <- measure_deriv / estimated_drymass_g #(X/h) / gDW -> X/(h gDW)
  
  dats$biolector$data <- addData(
    data = dats$biolector$data,
    ID = paste0(meas, "_rate"),
    dat = measure_deriv
  )
  
  # Add the new units and descriptors
  pos = length(dats$biolector$biolector_parameters$units) + 1
  
  dats$biolector$biolector_parameters$units = c(
    dats$biolector$biolector_parameters$units,
    sprintf("%s~gDW^-1", dats$biolector$biolector_parameters$units[meas]),
    sprintf("%s~h^-1~gDW^-1", dats$biolector$biolector_parameters$units[meas])
  )
  
  dats$biolector$biolector_parameters$descriptors = c(
    dats$biolector$biolector_parameters$descriptors,
    sprintf("norm.~%s", dats$biolector$biolector_parameters$descriptors[meas]),
    sprintf("norm.~%s~rate", dats$biolector$biolector_parameters$descriptors[meas])
  )
  
  # Adjust the names of the new units and descriptors
  names(dats$biolector$biolector_parameters$units)[c(pos, pos+1)] <-  
    names(dats$biolector$biolector_parameters$descriptors)[c(pos, pos+1)] <- 
    c(paste0(meas, "_norm"), paste0(meas, "_rate"))
  
  # Note that there are new measurements to be plotted
  measures = c(
    measures,
    paste0(meas, "_norm"),
    paste0(meas, "_rate")
  )
}

# Note that drymass and O2 consumption also have to be plotted
measures <- c(measures, c("drymass", "drymass_g", "O2_consumption"))

# Calculate the theoretical rate of O2 uptake into the medium and normalize it
O2_add_rate <- (O2PercToConc(100) - O2_conc) * dats$biolector$biolector_parameters$kLa # mmol/(l h)
O2_add_rate <- O2_add_rate / estimated_drymass_g # mmol/(l h gDW)

# Calculate the O2 consumption rate [mmol/(h gDW)]
# Here, as the difference between theoretical and observed alteration rate of O2 concentration
O2_consumption_rate <- 
  (getData(dats$biolector$data, "O2_conc_rate") - O2_add_rate) # mmol/(l h gDW)

O2_consumption_rate <- 
  - O2_consumption_rate * (dats$biolector$biolector_parameters$well_volume * 1e-3) # mmol/(h gDW), define consumption rate as positive

# Save the O2 consumption rate
dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "O2_consumption",
  dat = O2_consumption_rate
)

# Calculate the growthrates
for (meas in c("drymass", "scatter", "ribof")) {
  
  # Get the data to calculate growthrates
  measure_data <- getData(dats$biolector$data, meas)
  
  # Smooth the data with sm.spline
  measure_smooth <- smoothDat(
    varDat = measure_data,
    time = times,
    nOrder = 3,
    na.rm = T
  )$smthDat
  
  # Calculate the growthrates and save them
  temp_data <- addData(
    data = dats$biolector$data,
    ID = paste0(meas, "_smooth"),
    dat = measure_smooth
  )
  
  dpseg_model <- dpseg_plate(
    data = temp_data,
    yid = paste0(meas, "_smooth"), 
    P = .0001
  )
  
  dats$biolector$data <- addModel(
    fit = dpseg_model,
    data = dats$biolector$data,
    ID = sprintf("%s.mu.dpseg", meas),
    add.slopes = TRUE,
    col = "#FF0000"
  )
  
  # Add the new unit and descriptor
  pos <- length(dats$biolector$biolector_parameters$units)+1
  
  dats$biolector$biolector_parameters$units <- c(
    dats$biolector$biolector_parameters$units,
    "h^-1"
  )
  
  dats$biolector$biolector_parameters$descriptors <- c(
    dats$biolector$biolector_parameters$descriptors,
    sprintf("%s~growth~rate", dats$biolector$biolector_parameters$descriptors[meas])
  )
  
  # Adjust the name of the new unit and descriptor
  names(dats$biolector$biolector_parameters$units)[pos] <- 
    names(dats$biolector$biolector_parameters$descriptors)[pos] <- 
    sprintf("%s.mu.dpseg",meas)
  
  # Note that is a new measurement to be plotted
  measures <- c(measures,  sprintf("%s.mu.dpseg",meas))
}


# Plots of BioLector data  ---------------------------

# Set background colors
bgCols <- dats$biolector$bgCols

bgCols_df <- data.frame(
  col = bgCols,
  groups = factor_amount_uniq
)

# Create a label for the varied factor
factor_label <- sprintf("%s~'['*%s*']'", dats$biolector$varied_factor, dats$biolector$varied_factor_unit)
factor_label <- parse(text = factor_label)

# Save the background colors as a ggplot object
plts$objects[["bgCols"]] <- geom_rect(
  data = bgCols_df,
  aes(fill = col),
  fill = bgCols,
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.25
)

# Save the varied factor label as a ggplot object
plts$objects[["label_factor"]] <- list(
  scale_y_continuous(sec.axis = dup_axis(name = factor_label, labels = NULL)),
  theme(axis.ticks.y.right = element_blank())
)

# Save the point of interest annotations as a ggplot object
plts$objects[["poi_annotation"]] <- create_poi_annotation(dats$biolector$anno)

# Plot each measurement
for(meas in measures){

  # Get the measurement data
  measure_data <- getData(dats$biolector$data, meas)
  
  # Convert into a ggplot usable format
  measure_df <- ggplotDf(measure_data, factor_amount, times)
  measure_df <- measure_df[complete.cases(measure_df),]
  
  # Create a suitable y label
  y_label <- parse(
    text = sprintf(
      "%s~'['*%s*']'",
      dats$biolector$biolector_parameters$descriptors[meas],
      dats$biolector$biolector_parameters$units[meas]
    )
  )
  
  # Plot the measurement
  measure_plot = ggplot(measure_df) +
    theme_bw() +
    plts$objects$bgCols +
    plts$objects$poi_annotation +
    geom_line(aes(
      x = times,
      y = value,
      group = as.factor(labels)
    )) +
    
    scale_x_continuous(expand = c(0,0))+
    labs(x = "time [h]", y = y_label)
  
  # If the varied factor has more than one level the plot should be faceted
  facet_bool <- length(factor_amount_uniq) > 1
  if(facet_bool){
    measure_plot <- measure_plot + facet_grid(groups~.) + plts$objects$label_factor
  }      
  
  
  plts[[meas]] <- measure_plot
}

# Save the plots as files
savePlts(
  x = plts[names(plts) != "objects"],
  type = "pdf",
  filePath = FIGUREFOLDER,
  prefix = paste(
    dats$biolector$study_id,
    dats$biolector$study_token,
    sep = "_"
  ),
  width = 10,
  height = 6
)

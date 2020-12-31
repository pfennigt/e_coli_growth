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
study_id <- "SDY002"
study_sub <- NA

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


# Integrate drymass measurements  ---------------------------

# Define assay path
assay_id <- "ASY002"
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

# Get weights before ("pre") and after ("post") biomass application
weight_type = sub("(.*)_[0-9]+", "\\1", names(dat$data))
weight_type_uniqe = unique(weight_type[-1])

# Get weight sample time points
sample_times <- dat$metadata$time_h

# Plot raw data
sample_weights = reshape2::melt(
  dat$data,
  id.vars = "sample_id"
)

sample_weights$variable <- sub("(.*)_[0-9]+", "\\1", sample_weights$variable)
sample_weights$variable <- factor(
  sample_weights$variable, 
  levels = unique(sample_weights$variable)
)

plts$filter_weights_raw <- 
  ggplot(sample_weights, aes(x=NA, y=value, color=variable)) + 
    theme_bw() +
    geom_boxplot() +
    geom_point(position=position_dodge(0.75)) +
    facet_grid(.~factor(sample_id)) +
    labs(x="sample ID", y="filter weights [mg]", color="measurement") +
    scale_x_discrete(position="top", breaks=NULL)

# Select data collected after 2 days of drying for further analysis
analysis_weights <- dat$data[weight_type %in% c("sample_id", "pre", "post2")]

# Add the metadata
analysis_weights <- merge.data.frame(
  dat$metadata[,c("sample_id", "time_h", "volume_ml")],
  analysis_weights,
  by = "sample_id"
)

# Reshape
analysis_weights <- reshape2::melt(
  analysis_weights, 
  id.vars = c("sample_id", "time_h", "volume_ml"),
  value.name = "weight_mg",
  variable.name = "weight_type"
)

analysis_weights <- separate(
  analysis_weights,
  col = "weight_type",
  into = c("weight_type", "replicate"),
  sep = "_"
)

# Get mean "pre" weight of repeated measurements
pre_means = rowMeans(dat$data[weight_type == "pre"])
names(pre_means) = as.character(dat$metadata$sample_id)

# Track the estimated standard deviations
weight_sds <-  data.frame(sample_id = dat$data$sample_id)

weight_sds$pre <- apply(dat$data[weight_type == "pre"], 1, sd)

weight_sds$post2 <- apply(dat$data[weight_type == "post2"], 1, sd)

# Calculate weight differences
analysis_weights$diff_mg <- analysis_weights$weight_mg - pre_means

weight_sds$diff_mg <- sqrt(weight_sds$pre^2 + weight_sds$post2^2)

# Normalize weight differences
analysis_weights$norm_diff_mgml <- analysis_weights$diff_mg / analysis_weights$volume_ml

weight_sds$norm_diff_mgml <- weight_sds$diff_mg / dat$metadata$volume_ml

# Convert weight differences into carbon concentration
analysis_weights$drymass_cmmoll <- analysis_weights$norm_diff_mgml * cpg * mpc * 1000 #C-mmol/l

weight_sds$drymass_cmmoll <- weight_sds$norm_diff_mgml * cpg * mpc * 1000

# Plot the calculated dry mass
plts$drymass_weights <- 
  ggplot(analysis_weights[grepl("post", analysis_weights$weight_type),]) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype  = "dashed") +
    
    geom_point(aes(x = time_h, y = drymass_cmmoll)) +
    labs(x = "time [h]", y = "dry mass [C-mmol"~l^-1*"]", color = "Dryweight\nmeasurement")

# Get sampling times of wells
well_sampling_times = rep(Inf, length(wells))
names(well_sampling_times) <- wells

for (i in 1:nrow(dat$metadata)) {
  sampled_wells <- dat$metadata[i, "wells"]
  sampled_wells <- strsplit(sampled_wells,";")[[1]]
  
  well_sampling_times[sampled_wells] <- dat$metadata[i, "time_h"]
}

# Get BioLector measurements closest to the sampling times
nearest_bl_times <- data.table(time = times, val = times, ID = 1:length(times))
setattr(nearest_bl_times, "sorted", "time")  # let data.table know that w is sorted

nearest_bl_times <- nearest_bl_times[J(dat$metadata$time_h), roll = "nearest"]

nearest_bl_values <- lapply(c("scatter", "ribof"), function(x){
  getData(dats$biolector$data, x)[nearest_bl_times$ID, reference_wells]
})
names(nearest_bl_values) <- c("scatter", "ribof")

# Summarize data
summary_weights <- data.table(analysis_weights)
summary_weights <- summary_weights[
  weight_type == "post2",
  .(drymass_cmmoll = mean(drymass_cmmoll), norm_diff_mgml = mean(norm_diff_mgml)),
  by="time_h"
]
summary_weights <- as.data.frame(summary_weights)
summary_weights$sd_drymass_cmmoll <- weight_sds$drymass_cmmoll
summary_weights$sd_norm_diff_mgml <- weight_sds$norm_diff_mgml

summary_weights$scatter <- rowMeans(nearest_bl_values$scatter)
summary_weights$sd_scatter <- apply(nearest_bl_values$scatter, 1, sd)

summary_weights$ribof <- rowMeans(nearest_bl_values$ribof)
summary_weights$sd_ribof <- apply(nearest_bl_values$ribof, 1, sd)

# Annotate the approximate growth phase at the sampling points
summary_weights$growth_phase = c(rep("1. growth",4), rep("2. growth",2), "stationary")

# Plot the summarized data
summary_weights_plot = data.frame(
  time_h = summary_weights$time_h,
  norm_diff_mgml = summary_weights$norm_diff_mgml,
  sd_norm_diff_mgml = summary_weights$sd_norm_diff_mgml,
  drymass_cmmoll = summary_weights$drymass_cmmoll,
  sd_drymass_cmmoll = summary_weights$sd_drymass_cmmoll,
  growth_phase = summary_weights$growth_phase,
  
  value = c(summary_weights$scatter, summary_weights$ribof),
  sd_value = c(summary_weights$sd_scatter, summary_weights$sd_ribof),
  value_type = rep(c("scatter","ribof"), each = 7)
)

plts$drymass_vs_scatter_ribof <- 
  ggplot(summary_weights_plot) + 
    theme_bw() +
    geom_rect(
      data = data.frame(1),
      fill = dats$biolector$bgCols,
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.2
    ) +
    geom_rect(
      data = scatter_drop_values,
      aes(xmin = xmin, xmax = xmax, fill = factor(drop_no)),
      min = -Inf,
      ymax = Inf,
      alpha = 0.2
    ) +
    geom_hline(yintercept = 0, linetype ="dashed")+
    
    geom_point(aes(x = value, y = norm_diff_mgml, shape = growth_phase), size = 2.5)+
    geom_errorbar(
      aes(
        x = value,
        y = norm_diff_mgml,
        ymin = norm_diff_mgml - sd_norm_diff_mgml,
        ymax = norm_diff_mgml + sd_norm_diff_mgml
      ),
      width = 0.01) +
    geom_errorbar(
      aes(
        x = value,
        y = norm_diff_mgml,
        xmin = value - sd_value,
        xmax = value + sd_value
      ),
      width = 0.01
    ) + 
    
    facet_grid(
      value_type~.,
      labeller = labeller(value_type = c("ribof" = "riboflavin", "scatter" = "scatter"))
    ) +
    expand_limits(x = 0, y = 0)+
    labs(
      y = "dry mass [mg"~ml^-1*"]",
      x = "scatter or riboflavin [AU]",
      shape = "growth phase",
      fill = "scatter-drop"
    ) +
    scale_shape_manual(values = c(19, 15, 4))+
    scale_fill_manual(values = c("gray50", "mediumorchid4"), labels = list("No. 1", "No. 2")) +
    guides(
      shape = guide_legend(order = 1),
      fill = guide_legend(order = 2),
      color = guide_legend(order = 3)
    )

# Save the measured drymasses to add to the biomass plot later
plts$objects[["drymasses"]] = list(
  geom_point(
    data = summary_weights_plot,
    aes(
      x = time_h,
      y = norm_diff_mgml,
      shape = growth_phase
    ),
    size = 2.5,
    col = "mediumorchid4"),
  geom_errorbar(
    data = summary_weights_plot,
    aes(
      x = time_h,
      y = norm_diff_mgml,
      ymin = norm_diff_mgml - sd_norm_diff_mgml,
      ymax = norm_diff_mgml + sd_norm_diff_mgml
    ),
    width = 0.01,
    col = "mediumorchid4")
)

# Estimate drymass from riboflavin
drymass_lms <- list(
  description = "estimate dry mass [mg ml^-1] from riboflavin [AU]; values before the scatter drop with lm1, values after the scatter drop with lm2",
  estimator = "ribof",
  lm1 = lm(
    norm_diff_mgml ~ ribof,
    summary_weights[summary_weights$growth_phase == "1. growth",]
  ),
  lm2 = lm(
    norm_diff_mgml ~ ribof,
    summary_weights[summary_weights$growth_phase != "1. growth",]
  )
)

# Add linear models to the dry mass vs. riboflavin plot
# First create the appropritate label
drymass_lms_coefs <- lapply(drymass_lms[-c(1,2)], coef_na)
drymass_lms_coefs <- do.call(rbind, drymass_lms_coefs)
drymass_lms_coefs <- as.data.frame(drymass_lms_coefs)
colnames(drymass_lms_coefs) <-  c("intercept", "slope")
drymass_lms_coefs$value_type <- "ribof"

drymass_lms_label <- sprintf("%d. fitted slope: %.3f", c(1,2), drymass_lms_coefs$slope)
drymass_lms_label <- list(bquote(atop(.(drymass_lms_label[1]), .(drymass_lms_label[2]))))

# predict a few intermediate values to show in the plot
riboflavin_means = rowMeans(getData(dats$biolector$data, "ribof")[, reference_wells])

drymass_lms_predicts = data.frame(
  pred = estimate_two_phases(
    values = riboflavin_means,
    times = times,
    switchTime = scatter_drops[1, "t1"],
    lm1 = drymass_lms$lm1,
    lm2 = drymass_lms$lm2
  ),
  value = riboflavin_means,
  value_type = "ribof"
)

# Add the linear models to the plot
plts$drymass_vs_scatter_ribof <- 
  plts$drymass_vs_scatter_ribof +
    geom_abline(
      data = drymass_lms_coefs,
      aes(
        intercept = intercept,
        slope = slope,
        color = "growth phase\nspecific models"
      ),
      linetype = "dotted"
    ) +
    geom_label(
      data = data.frame(value_type = "ribof"),
      label = drymass_lms_label,
      x = Inf,
      y = 0,
      hjust = "right",
      vjust = "bottom",
      color = "black",
      parse = T
    ) +
    geom_path(data = drymass_lms_predicts, aes(x = value, y = pred, color = "concatenated")) +
    scale_color_manual(name="linear models", values=c("red", "black"))

# Save important variables in a list
dat$dataParameters <- list(
  sample_times <- dat$metadata$time_h,
  well_times <- well_sampling_times
)

# Save drymass data for later use
dats$drymass <- dat

# Save the estimation model to a file
drymass_lms_filename = paste(
  dat$study_id,
  dat$study_token,
  "drymass_estimation_twophase_ribof.rds",
  sep = "_"
)

saveRDS(
  drymass_lms,
  file = file.path(
    ANALYSISFOLDER,
    "data_analysis_outputs",
    drymass_lms_filename
  )
)


# Continue analysis of BioLector data  ---------------------------

# Create a mask excluding measured biolector values after sampling of the well
sampling_exclusion_mat = matrix(
  1,
  nrow = length(times),
  ncol = length(wells)
)

for (i in 1:ncol(sampling_exclusion_mat)) {
  sampling_exclusion_mat[
    c(
      dats$biolector$dataParameters$times[-1],
      Inf
    ) > well_sampling_times[i],
    i
  ] <- NA
} 

# Calculate the biomass concentration [gDW/l] and save
estimated_drymass_g <- apply(
  getData(dats$biolector$data, "ribof"),
  2,
  function(x){
    estimate_two_phases(
      values = x,
      times = times,
      switchTime = scatter_drops[1, "t1"],
      lm1 = drymass_lms$lm1,
      lm2 = drymass_lms$lm2
    )
  }
)

dats$biolector$data <- addData(
  data = dats$biolector$data,
  ID = "drymass_g",
  dat = estimated_drymass_g
)

# Convert to C-mmol/l and save
estimated_drymass <- estimated_drymass_g * cpg * mpc * 1000 # C-mmol/l

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

  # Get the measurement data and mask the sampled wells
  measure_data <- getData(dats$biolector$data, meas)
  measure_data <- measure_data * sampling_exclusion_mat
  
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

# add the measured dry masses to the respective plot
plts$drymass_g <- plts$drymass_g + plts$objects$drymasses

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

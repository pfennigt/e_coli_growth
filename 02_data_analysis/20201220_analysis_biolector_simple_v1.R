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
library(utils)
library(stringr)

# Set working directory to the investigation folder
INVESTIGATIONPATH <- "C:/Users/tpeng/Desktop/E_coli_growth"

biomassEqivalent = "Phase"

settleTime = 2 #0.5

kLa = 230 #h^-1

source()



#======================================= analysis =======================================

dats = list()
plts = list("glc"=list("items"=list()),
            "ace"=list("items"=list()),
            "glcBM"=list("items"=list()),
            "aceBM"=list("items"=list()),
            "glcRe1"=list("items"=list()),
            "aceRe1"=list("items"=list()),
            "glcInoc"=list("items"=list()),
            "glcSml"=list("items"=list()),
            "glcMdl" = list(),
            "all"=list())
dpsegs = list("glc"=list(),
              "ace"=list(),
              "glcRe1"=list(),
              "aceRe1"=list(),
              "glcInoc"=list(),
              "glcSml"=list())
coefs = list("glc"=list(),
             "ace"=list(),
             "glcRe1"=list(),
             "aceRe1"=list(),
             "glcInoc"=list(),
             "glcSml"=list(),
             "all"=list())
smooths = list("glc"=list(),
               "ace"=list(),
               "glcBM"=list(),
               "aceBM"=list(),
               "glcRe1"=list(),
               "aceRe1"=list(),
               "glcInoc"=list(),
               "glcSml"=list())
units = list("glc"=list(),
             "ace"=list(),
             "glcBM"=list(),
             "aceBM"=list(),
             "glcRe1"=list(),
             "aceRe1"=list(),
             "glcInoc"=list(),
             "glcSml"=list())
BMConvs = list("glcBM"=list(),
               "aceBM"=list())

#=================================== get biomass data ===================================
if(PC){INPATHHeader = "D://"}else{INPATHHeader = "C://Users/"}
source(paste0(INPATHHeader,'tpeng/OneDrive/OD__Universität/BA__Axmann/R/biolector_eval_biomass.R'), echo=TRUE)



#=========================================================================================

INPATH = paste0(INPATHHeader,INPATH)

DATPATH <- file.path(INPATH,"data")
RESPATH <- file.path(INPATH,"results") # results

setwd(RESPATH)


# expid <- "25_Ecoli_2020_REDUCTION-1"
# data.file <- file.path(DATPATH, paste0(expid,".csv"))
# layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
# 
# 
# #view full dataset
# dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
#                       layout=layout.file, 
#                       fields=c("strain","glc","ace","aTc"), 
#                       afields=c("glc","ace","aTc"),blank.id="blank",
#                       group1="glc", group2 = c("glc","amount"),
#                       group2.color="color")

# viewPlate(dat)
# 
# viewMap(map=dat$layout, nrow=6, ncol=8, color="color", text="amount", title="plate layout: glucose")
# viewMap(map=dat$layout, nrow=6, ncol=8, color="ace.color", text="ace.amount", title="plate layout: acetate")



# CsourceShort = "glc"
# CsourceShort = "ace"
# CsourceShort = "glcRe1"
# CsourceShort = "aceRe1"
# CsourceShort = "glcSml"
# CsourceShort = "glcInoc"

for(CsourceShort in c("glc", "ace", "glcRe1", "aceRe1", "glcSml", "glcInoc")){
  
  if(CsourceShort =="glc"){
    CSource = "glucose"
    CSourceRead = "glc"
    
    expid <- "25_Ecoli_2020_REDUCTION-1"
    data.file <- file.path(DATPATH, paste0(expid,".csv"))
    layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
    
    skipWellList = c(paste0(toupper(letters[4:6]), rep(1:7,3)))
    Amounts = "amount"
    PltColors = "color"
    bgCols = colorRampPalette(c("white", "royalblue3"))(7)
    xrng=c(0,12.5)
    
    ignoreWells = paste0(rep(LETTERS[1:6], each=2), c(1,8))
    peakTimeRange = c(5.5, 12.5)
    
    
    
    #read data
    dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                          layout=layout.file, 
                          fields=c("strain","glc","ace","aTc"), 
                          afields=c("glc","ace","aTc"),
                          blank.id="blank",blank.data=c("Biomass","Riboflavine", "NADH - NADPH"),
                          skip.wells = skipWellList,
                          group1=CSourceRead, group2 = c(CSourceRead,Amounts),
                          group2.color=PltColors)
    
  }else if(CsourceShort =="ace"){
    CSource = "acetate"
    CSourceRead = "ace"
    
    expid <- "25_Ecoli_2020_REDUCTION-1"
    data.file <- file.path(DATPATH, paste0(expid,".csv"))
    layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
    
    skipWellList = c(paste0(toupper(letters[1:3]), rep(1:7,3))#, "E7" #additionally E7, sinceit shows too large scatter values
                     ) #additionally E7, sinceit shows too large scatter values
    Amounts = "ace.amount"
    PltColors = "ace.color"
    bgCols = colorRampPalette(c("white", "firebrick2"))(7)
    xrng=c(0,37.5)
    
    ignoreWells = paste0(rep(LETTERS[1:6], each=2), c(1,8))
    peakTimeRange = c(15, 25)
    
    #scatterPeakThresh = -0.02
    #scatterPeakMinL = 5 # slopes are long and gentle
    #scatterPeakThresh = -0.015
    #scatterPeakMinL = 15 # slopes are long and gentle
    
    #read data
    dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                          layout=layout.file, 
                          fields=c("strain","glc","ace","aTc"), 
                          afields=c("glc","ace","aTc"),
                          blank.id="blank",blank.data=c("Biomass","Riboflavine", "NADH - NADPH"),
                          skip.wells = skipWellList,
                          group1=CSourceRead, group2 = c(CSourceRead,Amounts),
                          group2.color=PltColors)
    
  }else if(CsourceShort == "glcRe1"){
    CSource = "glucose"
    CSourceRead = "glc"
    
    expid = "28_Ecoli_2020_REDUCTION-1"
    data.file <- file.path(DATPATH, paste0(expid,".csv"))
    layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
    
    Amounts = "amount"
    PltColors = "color"
    bgCols = colorRampPalette(c("white", "royalblue3"))(7)
    xrng=c(0,12.5)
    
    ignoreWells = c(paste0(LETTERS[1:6], 1), "D8","E8","F8")
    peakTimeRange = c(5.5, 12.5)
    
    #scatterPeakThresh = -0.1
    #scatterPeakMinL = 5
    
    #read data
    dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                          layout=layout.file, 
                          fields=c("strain","glc","ace","aTc"), 
                          afields=c("glc","ace","aTc"),
                          blank.id="blank",blank.data=c("Biomass","Riboflavine", "NADH - NADPH"),
                          group1=CSourceRead, group2 = c(CSourceRead,Amounts),
                          group2.color=PltColors)
    
  }else if(CsourceShort == "aceRe1"){
    CSource = "acetate"
    CSourceRead = "ace"
    
    expid = "34_Ecoli_2020-1831_REDUCTION-1"
    data.file <- file.path(DATPATH, paste0(expid,".csv"))
    layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
    
    Amounts = "ace.amount"
    PltColors = "ace.color"
    bgCols = colorRampPalette(c("white", "firebrick2"))(7)
    xrng=c(0,60)
    
    ignoreWells = c(paste0(LETTERS[1:6], 1), "D8","E8","F8")
    peakTimeRange = c(30,50)
    
    #scatterPeakThresh = -0.1
    #scatterPeakMinL = 5
    
    #read data
    dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                          layout=layout.file, 
                          fields=c("strain","glc","ace","aTc"), 
                          afields=c("glc","ace","aTc"),
                          blank.id="blank",blank.data=c("Biomass","Riboflavine", "NADH - NADPH"),
                          group1=CSourceRead, group2 = c(CSourceRead,Amounts),
                          group2.color=PltColors)
    
  }else if(CsourceShort == "glcSml"){
    CSource = "glucose"
    CSourceRead = "glc"
    
    expid = "32_Ecoli_2020-1831_REDUCTION-1"
    data.file <- file.path(DATPATH, paste0(expid,".csv"))
    layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
    
    Amounts = "amount"
    PltColors = "color"
    bgCols = colorRampPalette(c("white", "royalblue3"))(15)
    xrng=c(0,15)
    
    ignoreWells = c(paste0(rep(LETTERS[1:3], each=4), 1:4), "D8","E8","F8")
    peakTimeRange = c(5.5, 14.5)
    
    #scatterPeakThresh = -0.1
    #scatterPeakMinL = 5
    
    #read data
    dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                          layout=layout.file, 
                          fields=c("strain","glc","ace","aTc"), 
                          afields=c("glc","ace","aTc"),
                          blank.id="blank",blank.data=c("Biomass","Riboflavine", "NADH - NADPH"),
                          group1=CSourceRead, group2 = c(CSourceRead,Amounts),
                          group2.color=PltColors)
    
    
  }else if(CsourceShort == "glcInoc"){
    CSource = "glucose"
    CSourceRead = "glc"
    
    expid = "33_Ecoli_2020-1831_REDUCTION-1"
    data.file <- file.path(DATPATH, paste0(expid,".csv"))
    layout.file <- file.path(DATPATH, paste0(expid,"_layout.csv"))
    
    #skipWellList = "C1" #C1 shows abnormal initial ascent
    Amounts = "inoc.amount"
    PltColors = "inoc.color"
    bgCols = colorRampPalette(c("royalblue3", "green4"))(500)[c(20, 25, 32, 40, 50, 63, 79, 100, 126, 158, 199, 251, 316, 397, 500)]
    xrng=c(0,18.5)
    
    ignoreWells = c("D8","E8","F8")
    peakTimeRange = c(5, 17)
    
    #scatterPeakThresh = -0.1
    #scatterPeakMinL = 5
    
    #read data
    dat <- readExperiment(data.file, type="BioLectorPro", time.conversion=1/3600,
                          layout=layout.file, 
                          fields=c("strain","glc","ace","aTc","inoc"), 
                          afields=c("glc","ace","aTc","inoc"),
                          blank.id="blank",blank.data=c("Biomass","Riboflavine", "NADH - NADPH"),
                          #skip.wells = skipWellList,
                          group1=CSourceRead, group2 = c(CSourceRead,Amounts),
                          group2.color=PltColors)
  }
  
  dat <- prettyData(dat, yids=c(ribof="Riboflavine",O2="DO(Pst3)",
                                scatter="Biomass", pH="pH(HP8)",
                                NADH="NADH - NADPH"))
  
  #format data
  #dat$layout[[Amounts]] <- (dat$layout[[Amounts]]/max(dat$layout[[Amounts]])) * MaxConc[CsourceShort] #convert to C-mmol
  #layout already in C-mM
  dat$layout[[Amounts]][dat$layout$strain=="blank"] <- NA
  dat$layout$color[dat$layout$strain=="blank"] <- "#999999"
  
  
  #cut off long stationary phase
  # xStat = findSatGrouped(dat, "scatter", Amounts, npoints = 10)
  
  dat <- cutData(dat, xrng=xrng)
  nTime = length(dat$Time)

  concs = dat$layout[[Amounts]][order(dat$layout$well)]
  concsUniq = sort(unique(concs))
  times = dat$Time
  wells = sort(dat$layout$well)
  concCols = c("#BEBEBE", bgCols[-1]) #for coloring points: 0 C-mmol is grey, not white
  bgColsDf = data.frame(col=bgCols, groups = sort(unique(dat$layout[,Amounts])))
  
  xLabel = if(CsourceShort == "glcInoc"){
    "inoculation OD"
  }else{
    sprintf("%s concentration [C-mM]", CSource)
  }
  
  dat$concs = concs
  
  #adjust O2
  O2 = getData(dat, "O2")
  O2 = O2/max(O2)*100
  
  dat = addData(dat, "O2", O2, replace = T)
  
  
  #POI anootation data
  POIDf = read.csv(file.path(INPATH,"data", "POI", paste0(CsourceShort,"_POI.csv")), header = T)
  POIDf$group = as.character(POIDf$group)
  POIDf$pos = as.numeric(paste0(POIDf$pos,"Inf"))
  POIDf$vjust = mapvalues(x = POIDf$pos, to=c("bottom","top"), from = c(-Inf, Inf))
  
  #raw plots
  # viewGroups(dat, yids="scatter", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["scatter"]]= recordPlot()
  # viewGroups(dat, yids="ribof", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["ribof"]]= recordPlot()
  # viewGroups(dat, yids="O2", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["O2"]]= recordPlot()
  # viewGroups(dat, yids="pH", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["pH"]]= recordPlot()
  # viewGroups(dat, yids="NADH", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["NADH"]]= recordPlot()
  
  
  
  
  #ggplot mods
  bgColsGG = geom_rect(data = bgColsDf, aes(fill = col), fill= bgCols, xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.25)
  
  timeStartGG = geom_vline(xintercept = settleTime, linetype="dotted")
  

  labelFacetGG = list(scale_y_continuous(sec.axis = dup_axis(name = xLabel, labels=NULL)),
                        theme(axis.ticks.y.right = element_blank()))

  
  POIanno = list(geom_point(data = POIDf, aes(x=x, y=pos), size=8,shape=21, fill="grey"),
                 geom_text(data=POIDf, aes(x=x, y=pos, vjust=vjust, label=lab),size=3.1))
  
  
  
  
  # #find peaks via dpseg and mark them in plots
  # scatterMean = sapply(concsUniq, function(x){
  #   rowMeans(getData(dat, "scatter")[,concs == x], na.rm = T)
  # })
  # 
  # #smoothing necessary?
  # scatterMean = smoothDat(scatterMean, dat$Time)$smthDat
  # # scatterMean = apply(scatterMean, 2, ma, 3)
  # 
  # scatterPeaksDf = findDrops(mat = scatterMean,
  #                            times = dat$Time,
  #                            groups = concsUniq,
  #                            slopeThresh = scatterPeakThresh,
  #                            minl = scatterPeakMinL,
  #                            timeThresh = 1)
  # 
  # #more
  # 
  # scatterPeaksGG = geom_rect(data = scatterPeaksDf, aes(xmin = t1, xmax = t2), fill= "grey", ymin = -Inf,ymax = Inf, alpha = 0.4)
  
  
  
  #find drops and their times via <5% normal - slopes
  scatterPeaksAllDf = 
    findDrops(mat = getData(dat, "scatter"), #find drops in scatter
              times = dat$Time,
              groups = concs,
              derivQuantThresh = 0.05,
              derivQuantThreshExpand = 0.1,
              inclTimeRange = settleTime,
              exclLast = T)
  
  scatterPeaksAllDf$well = dat$layout$well[order(dat$layout$well)][scatterPeaksAllDf$column]
  scatterPeaksAllDf$row = sub("([A-Z]+)([1-9])", "\\1", scatterPeaksAllDf$well)
  scatterPeaksAllDf$col = sub("([A-Z]+)([1-9])", "\\2", scatterPeaksAllDf$well)
  
  #show peaks in wells
  scatter = getData(dat,"scatter")
  
  scatterPlateDf = ggplotDf(scatter, concs, times)
  #scatterPlateDf = scatterPlateDf[complete.cases(scatterPlateDf),]
  scatterPlateDf$row = sub("([A-Z]+)([1-9])", "\\1", scatterPlateDf$labels)
  scatterPlateDf$col = sub("([A-Z]+)([1-9])", "\\2", scatterPlateDf$labels)
  
  if(CsourceShort %in% c("glc", "ace")){
      bgColsPlate = data.frame(row=if(CsourceShort == "glc"){c("A","B","C")}else if(CsourceShort == "ace"){c("D","E","F")},
                           col=rep(1:7, each=3),
                           cols=rep(bgCols, each=3))
      
  }else if(CsourceShort %in% c("glcRe1", "aceRe1")){
    bgColsPlate = data.frame(row=LETTERS[1:6],
                             col=rep(1:7, each=6),
                             cols=rep(bgCols, each=6))
    bgColsPlate = rbind(bgColsPlate, data.frame(row=LETTERS[1:3], col=8, cols=bgCols[7]))
    
  }else if(CsourceShort =="glcSml"){
    bgColsPlate = data.frame(row=rep(c("A","D","B","E","C","F"), each=8),
                             col=1:8,
                             cols=rep(bgCols[c(1:15,NA)], 3))
    bgColsPlate[is.na(bgColsPlate)] = bgCols[c(1,8,15)]
    
  }else if(CsourceShort =="glcInoc"){
    bgColsPlate = data.frame(row=rep(c("A","D","B","E","C","F"), each=8),
                             col=1:8,
                             cols=rep(c(bgCols[1:15],"royalblue3"), 3))
  }

  
  #scatterPlateDf = scatterPlateDf[complete.cases(scatterPlateDf),]
  scatterPlatePlt = ggplot(scatterPlateDf)+theme_bw()+
    geom_rect(data = bgColsPlate, aes(fill = cols), fill= bgColsPlate$cols, xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf, alpha = 0.25)+
    geom_rect(data = scatterPeaksAllDf, aes(xmin = t1, xmax = t2), fill= "grey30", ymin = -Inf,ymax = Inf, alpha = 0.4)+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(row~col)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y="scatter [AU]")
  plts[[CsourceShort]][["scatterPeaksPlateGG"]] = scatterPlatePlt
  
  
  #exclude artifacts outside of the expected range
  scatterPeaksAllDf = scatterPeaksAllDf[!(scatterPeaksAllDf$t2<peakTimeRange[1] | scatterPeaksAllDf$t1>peakTimeRange[2]),]
  
  
  #summarize Peaks for groups by taking largest ranges for all peaks
  #dont use no-Carbon or blank columns
  scatterPeaksAllDf = scatterPeaksAllDf[!scatterPeaksAllDf$well %in% ignoreWells,]
  scatterPeaksAllDf_appear = table(scatterPeaksAllDf$column)
  
  scatterPeaksDf = lapply(unique(scatterPeaksAllDf$groups), function(conc){
    submat = scatterPeaksAllDf[scatterPeaksAllDf$groups==conc, c("t1", "t2", "n1", "n2", "groups")]
    submat = submat[order(submat$n1),]
    
    if(nrow(submat)<=1) return(submat)
    
    peakMat = submat[1,]
    
    for(sRow in 2:(nrow(submat))){
      subRow = submat[sRow,]
      brokenBool = F
      
      for(pRow in 1:nrow(peakMat)){ #do the ranges overlap?
        underBool = all(peakMat[pRow, 1:2] < subRow[,1])
        overBool = all(peakMat[pRow, 1:2] > subRow[,2])
        
        if(!(underBool | overBool)){
          peakMat[pRow, 1] = min(peakMat[pRow, 1], subRow[,1])
          peakMat[pRow, 3] = min(peakMat[pRow, 3], subRow[,3])
          peakMat[pRow, 2] = max(peakMat[pRow, 2], subRow[,2])
          peakMat[pRow, 4] = max(peakMat[pRow, 4], subRow[,4])
          
          brokenBool = T
          break
        }
      }
      if(!brokenBool){
        peakMat = rbind(peakMat, subRow)
      }
    }
    return(peakMat)
  })
  
  scatterPeaksDf = do.call(rbind, scatterPeaksDf)
  
  dat$scatterPeaks = scatterPeaksDf
  
  #mark peaks in following plots
  scatterPeaksGG = geom_rect(data = scatterPeaksDf, aes(xmin = t1, xmax = t2), fill= "grey", ymin = -Inf,ymax = Inf, alpha = 0.5)
  
  
  
  #Plos as boxplots for comparing peak times; only take the peaks beginning times for boxplots
  scatterPeaksPos = matrix(NA,length(scatterPeaksAllDf_appear),max(scatterPeaksAllDf_appear))
  row.names(scatterPeaksPos) = names(scatterPeaksAllDf_appear)
  
  for(k in names(scatterPeaksAllDf_appear)){
    columnPos = which(scatterPeaksAllDf$column == k)
    scatterPeaksPos[k,1:length(columnPos)] = scatterPeaksAllDf[columnPos, "t1"]
  }
  
  scatterPeaksPos = data.frame(scatterPeaksPos)
  colnames(scatterPeaksPos) = sub("X", "Peak",colnames(scatterPeaksPos))
  
  #if there are two or more peaks include peak time difference
  if(ncol(scatterPeaksPos)>=2){
    #also mark differences in peak times per well
    scatterPeaksPos_difs = sapply(2:ncol(scatterPeaksPos), function(colNum){
      scatterPeaksPos[,colNum] - scatterPeaksPos[,colNum-1]
    })
    
    #colnames(scatterPeaksPos_difs) = sapply(2:ncol(scatterPeaksPos), function(colNum){paste(colnames(scatterPeaksPos)[c(colNum-1,colNum)], collapse="_")})
    colnames(scatterPeaksPos_difs) = paste0("Delta*", colnames(scatterPeaksPos)[-1])
    
    scatterPeaksPos = cbind(scatterPeaksPos, scatterPeaksPos_difs)
  }
  
  
  #add the respectibe glucose level
  scatterPeaksPos$groups = sapply(names(scatterPeaksAllDf_appear), function(x){
    scatterPeaksAllDf[scatterPeaksAllDf$column == x, "groups"][1]
  })
  
  scatterPeaksPosDf = reshape2::melt(scatterPeaksPos, id.vars = "groups")
  scatterPeaksPosDf = scatterPeaksPosDf[complete.cases(scatterPeaksPosDf),]
  
  
  bgColsBox = (concsUniq[-1] - concsUniq[-length(concsUniq)])/2
  bgColsBox1 = bgColsBox[c(1, 1:length(bgColsBox))]
  bgColsBox2 = bgColsBox[c(1:length(bgColsBox), length(bgColsBox))]
  bgColsBoxDf = data.frame(xmin=concsUniq - bgColsBox1, xmax=concsUniq+bgColsBox2, cols= bgColsGG$data$col)
  
  scatterPeaksPosPlt = ggplot()+theme_bw()+
    geom_rect(data = bgColsBoxDf, aes(xmin =xmin, xmax=xmax, fill = cols), fill= rep(bgColsBoxDf$cols, length(unique(scatterPeaksPosDf$variable))),ymin = -Inf,ymax = Inf, alpha = 0.25)+
    geom_hline(yintercept = 0, linetype ="dashed")+
    
    geom_point(data=scatterPeaksPosDf, aes(x=groups, y=value, group=groups))+
    
    labs(y="time of respective drops beginning [h]", x=xLabel)+
    scale_x_continuous(expand=c(0,0))
  
  plts[[CsourceShort]][["scatterPeaksPosGG"]] = scatterPeaksPosPlt + if(length(unique(scatterPeaksPosDf$variable))>=2) {facet_grid(variable~., labeller = label_parsed)}
  
  
  
  #times of manually annotated points of interest
  POIPosPltDf = POIDf
  if(CSourceRead == "glc"){
    POILabs = c("EOL", "FHP", "SD1", "SD2")
    POIPosPltDf$lab = factor(mapvalues(POIPosPltDf$lab, to=POILabs, from = as.character(1:4)), levels = POILabs)
  }else{
    POILabs=c("EOL", "GRI", "O2D", "GRD", "SDA")
    POIPosPltDf$lab = factor(mapvalues(POIPosPltDf$lab, to=POILabs, from = as.character(1:5)), levels = POILabs)
  }
  
  
  POIPosPlt = ggplot()+theme_bw()+
    geom_rect(data = bgColsBoxDf, aes(xmin =xmin, xmax=xmax, fill = cols), fill= rep(bgColsBoxDf$cols, length(unique(POIDf$lab))),ymin = -Inf,ymax = Inf, alpha = 0.25)+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    
    geom_point(data=POIPosPltDf, aes(x=groups, y=x, group=groups))+
    
    labs(y="time of respective POI [h]", x=xLabel)+
    scale_x_continuous(expand=c(0,0))
  plts[[CsourceShort]][["POIPosGG"]] = POIPosPlt + if(length(unique(POIDf$lab))>=2) {facet_grid(lab~., scales = "free_y")}
  
  
  #times of manually annotated phases
  POIList = split.data.frame(POIDf, f = POIDf$lab)
  POIList[["0"]] = POIList[["1"]]
  POIList[["0"]]$x = 0
  
  if(CSourceRead == "glc"){
    PhaseLabs = c("lag phase", "slow first\ngrowth phase", "fast first\ngrowth phase ", "second\ngrowth phase ")
    PhaseList = lapply(as.character(1:4), function(lab){
      res = POIList[[lab]]
      res$x = res$x - POIList[[as.character(as.numeric(lab)-1)]]$x
      return(res)
    })
    PhaseDf = do.call(rbind, PhaseList)
    PhaseDf$lab = factor(mapvalues(PhaseDf$lab, to=PhaseLabs, from = as.character(1:4)), levels = PhaseLabs)
    
  }else{
    PhaseLabs = c("lag phase", "slow \ngrowth phase", "fast \ngrowth phase", bquote(atop("reduced"~O[2],"growth phase")))
    PhaseList = lapply(as.character(1:4), function(lab){
      res = POIList[[lab]]
      res$x = res$x - POIList[[as.character(as.numeric(lab)-1)]]$x
      return(res)
    })
    PhaseDf = do.call(rbind, PhaseList)
    PhaseDf$lab = factor(mapvalues(PhaseDf$lab, to=PhaseLabs, from = as.character(1:4)), levels = PhaseLabs)
  }
  
  PhasePlt = ggplot()+theme_bw()+
    geom_rect(data = bgColsBoxDf, aes(xmin =xmin, xmax=xmax, fill = cols), fill= rep(bgColsBoxDf$cols, length(unique(PhaseDf$lab))),ymin = -Inf,ymax = Inf, alpha = 0.25)+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    
    geom_point(data=PhaseDf, aes(x=groups, y=x, group=groups))+
    
    labs(y="length of respective phase [h]", x=xLabel)+
    scale_x_continuous(expand=c(0,0))
  plts[[CsourceShort]][["PhaseLengthGG"]] = PhasePlt + if(length(unique(PhaseDf$lab))>=2) {facet_grid(lab~., scales = "free_y",labeller = label_parsedMixed)}
  
  
  
  
  
  plts[[CsourceShort]][["items"]][["stackPlotBasis"]]=
    ggplot()+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    
    facet_grid(groups~.)+
    labs(x="time [h]", y="mean percentage of absolute max [%]", color="measured")+
    scale_x_continuous(expand=c(0,0))
  
  #make ggplots for all measured variables and create smoothed counterparts
  for(j in list(list("scatter","scatter [AU]"),
                list("ribof","riboflavin [AU]"),
                list("O2",bquote(O[2]~"saturation [%]")),
                list("pH","pH"),
                list("NADH","NADH [AU]"))){
  
    varDat = getData(dat,j[[1]])
    
    pltDf = ggplotDf(varDat, concs, times)
    pltDf = pltDf[complete.cases(pltDf),]
    pltPlt = ggplot(pltDf)+theme_bw()+
      bgColsGG+
      timeStartGG+
      labelFacetGG+
      scatterPeaksGG+
      POIanno+
      
      geom_line(aes(x=times, y=value, group=as.factor(labels)))+
      
      facet_grid(groups~.)+
      scale_x_continuous(expand=c(0,0))+
      labs(x="time [h]", y=j[[2]])
    
    if(j[1] %in% c("scatter", "ribof", "NADH")) pltPlt = pltPlt + geom_hline(yintercept = 0, linetype ="dashed")
    
    plts[[CsourceShort]][[paste0(j[[1]],"GG")]] = pltPlt + if(CsourceShort == "ace" & j[1] == "pH"){coord_cartesian(ylim = c(5.7,7.5))}
    
    plts[[CsourceShort]][["items"]][[paste0(j[[1]],"GG")]] = geom_line(data = pltDf, aes(x=times, y=value, group=as.factor(labels)))
    
    
    #relative data for comparison
    pltDf = data.table(pltDf)
    pltDf = pltDf[,.(value=mean(value)), by=c("groups", "times")]
    pltDf = data2Perc(pltDf, settleTime)
    
    plts[[CsourceShort]][["items"]][[paste0(j[1],"StackGG")]] = geom_line(data = pltDf, aes_(x=quote(times), y=quote(value), group="mean", color=j[[1]]))
    
      
    #smooth variable
    varSM = smoothDat(varDat, dat$Time, nOrder=3)
    smooths[[CsourceShort]][[j[[1]]]] = varSM$smthObj #save smoothing data
    dat <- addData(dat, paste0(j[1],"SM"), dat=varSM$smthDat)
  }
  
  #save units for eventual lookup
  for( k in 1:12){
    kNam = c("scatter","ribof","O2","pH","NADH")[k]
    units[[CsourceShort]][[kNam]] = c("AU", "AU", "%", "", "AU")[k]
  }
  
  # #find scatter peaks
  # scatterPeak = getData(dat, "scatter")
  # scatterPeak = apply(scatterPeak, 2, ma, 5)
  # #is the following value bigger (1) or smaller (-1)
  # scatterPeak = (rbind(scatterPeak[-1,], NA) - scatterPeak) 
  # scatterPeak = sign(scatterPeak)
  # 
  # #changes in this sign show peaks or valleys
  # scatterPeak = (rbind(NA, scatterPeak[-nrow(scatterPeak),]) + scatterPeak)
  # row.names(scatterPeak) = names(dat$Time)
  
  # scatterSM = smooths[[CsourceShort]][["scatter"]]
  # scatterDerif = sapply(scatterSM, function(smthdat){
  #   predict(smthdat, xarg = dat$Time ,nderiv = 1)
  #   })
  # scatterDerif = sapply(concsUniq, function(x){
  #   rowMeans(scatterDerif[,concs == x], na.rm = T)
  # })
  # 
  # 
  # scatterDerif2 = sapply(scatterSM, function(smthdat){
  #   predict(smthdat, xarg = dat$Time ,nderiv = 2)
  # })
  # scatterDerif2 = sapply(concsUniq, function(x){
  #   rowMeans(scatterDerif2[,concs == x], na.rm = T)
  # })
  
  
  
  ribof = getData(dat, "ribof")
  scatter = getData(dat, "scatter") # a.u., arbitrary unit
  gps <- 0.603 # mg/ml/scatter, biomass per scatter
  cpg <- 0.474 # g/g, gram carbon per gram biomass, Folsom&Carlson 2015
  mpc <- 1/12 # mol/g, mol carbon per gram carbon
  
  
  #show realtion riboflavin ~ ~scatter
  RPerSDf = ggplotDf(ribof, concs, times)
  RPerSDf$value = RPerSDf$value/ ggplotDf(scatter, concs, times)$value
  RPerSDf = RPerSDf[complete.cases(RPerSDf),]
  RPerSPlt = ggplot(RPerSDf)+theme_bw()+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y="fraction riboflavin / scatter")
  plts[[CsourceShort]][["ribofPerScatGG"]]=RPerSPlt
  
  
  
  
  if(biomassEqivalent == "scatter"){
    #convert scatter to biomass
    biomass <- scatter * 1e3 * gps * cpg * mpc # mmol/L
    dat <- addData(dat, "biomass", dat=biomass)
  
    
  }else if(biomassEqivalent == "ribof"){
    #convert riboflavin to biomass
    
    #find stable and linear ("saturated") areas near the end of measurement for scatter and riboflavin for each group of C-molarity
    SatBeginnings = data.frame(
      ribof = findSatGrouped(dat, "ribof", Amounts),
      scatter = findSatGrouped(dat, "scatter", Amounts)
    )
    
    SatBeginnings = apply(SatBeginnings, 1, max) #take the later saturation for each group
    
    #Each group of Concentrations evaluated at smaller saturation range (ribof or scatter)
    #SatBeginnings = rep(SatBeginnings,each = 3)[order(dat$layout$well)]
    
    #Alternatively: All groups evaluated in smallest saturation range
    SatBeginnings = rep(max(SatBeginnings), 3*length(SatBeginnings))[order(dat$layout$well)]
    
    
    ribof_SatVec = lapply(1:length(SatBeginnings), function(x){
      if(is.na(SatBeginnings[x])){return(NULL)}
      
      ribof[SatBeginnings[x]:nTime,x]
    })
    ribof_SatVec=unlist(ribof_SatVec)
    
    scatter_SatVec = lapply(1:length(SatBeginnings), function(x){
      if(is.na(SatBeginnings[x])){return(NULL)}
      
      scatter[SatBeginnings[x]:nTime,x]
    })
    scatter_SatVec=unlist(scatter_SatVec)
    
    ribofScatterLm = lm(ribof_SatVec~0+scatter_SatVec)
    
    ribofScatterScale = coef(ribofScatterLm)
    ribofScatterRsq = summary(ribofScatterLm)$r.squared
    
    ribofScatterDf = data.frame(ribof=ribof_SatVec, scatter =scatter_SatVec)
    
    
    ribofScatterPlt_label = c(sprintf("slope: %.3f", ribofScatterScale), sprintf(": %.3f",ribofScatterRsq))
    ribofScatterPlt_label = list(bquote(atop(.(ribofScatterPlt_label[1]),~~~~R^2*.(ribofScatterPlt_label[2]))))
    
    
    ribofScatterPlt = ggplot()+theme_bw()+
      geom_point(data=ribofScatterDf, aes(x=scatter, y=ribof))+
      geom_point(aes(x=0,y=0), color="red", size=2)+
      geom_abline(slope=ribofScatterScale, color="red", linetype="dashed")+
      
      labs(x="scatter [AU]", y="riboflavin [AU]")+
      geom_label(label=ribofScatterPlt_label, aes(x=max(scatter_SatVec),y=0.5*max(ribof_SatVec)), hjust="right", color="red", parse=T)
    
    plts[[CsourceShort]][["ribof_scatterGG"]]= ribofScatterPlt
    
    biomass <- ribof/ ribofScatterScale * 1e3 * gps * cpg * mpc   # C-mmol/L, scale riboflavin to scater and convert to biomass with scatter scaling
    dat <- addData(dat, "biomass", dat=biomass)
  
  }else if(biomassEqivalent == "ribofLm"){
    biomass <- predictBM(getData(dat, "ribof"), BMConvs$ribofLm) * 1e3 *cpg * mpc   # C-mmol/l
    dat <- addData(dat, "biomass", dat=biomass)
  
  }else if(biomassEqivalent == "Phase"){
    #use growth phase dependent biomass models
    
    if(CSourceRead == "glc"){
      biomass = lapply(wells, function(well){
        RibofBMSwitchGroup = dat$layout[dat$layout$well == well, Amounts]
        
        if(is.na(RibofBMSwitchGroup)){
          RibofBMSwitchTime = Inf
        }else if(RibofBMSwitchGroup==0){
          RibofBMSwitchTime = Inf
        }else{
          RibofBMSwitchTime = scatterPeaksDf[scatterPeaksDf$groups == RibofBMSwitchGroup, "t1"][1]
          
          if(is.na(RibofBMSwitchTime)) RibofBMSwitchTime = Inf
        }


        RibofBM = BMConvs$glcBM$ribofPhaseLm$mdl(values = ribof[,well],
                                                 times,
                                                 switchTime = RibofBMSwitchTime,
                                                 lm1 = BMConvs$glcBM$ribofPhaseLm$lm1,
                                                 lm2 = BMConvs$glcBM$ribofPhaseLm$lm2
                                  )
        
        return(RibofBM)
      })
      biomass = do.call(cbind, biomass)
      colnames(biomass) = wells
      
    }else if(CSourceRead == "ace"){
      biomass = apply(scatter, 2, function(x){BMConvs$aceBM$scatterPhaseLm$mdl(values = x,
                                                                             times,
                                                                             switchTime = Inf,
                                                                             lm1 = BMConvs$aceBM$scatterPhaseLm$lm1,
                                                                             lm2 = BMConvs$aceBM$scatterPhaseLm$lm2
        )
      })
      
      #estimate only valid until first drop, replace everything after with NA
      cutOffTimes = sapply(concs, function(conc){
        tim = scatterPeaksDf[scatterPeaksDf$groups == conc, "t1"]
        if(length(tim >1)){
          tim = tim[1]
        }else if(length(tim)==0){
          tim = min(scatterPeaksDf$t1)
          return(tim)
        }
        
        if(is.na(tim)) tim = min(scatterPeaksDf$t1)
        return(tim)
      })
      
      cutOffBool = sapply(cutOffTimes, function(tim){times>tim})
      
      biomass[cutOffBool] = NA
    }
    
    biomass = biomass * 1e3 *cpg * mpc # C-mmol/l
    dat <- addData(dat, "biomass", dat=biomass)
  }
  
  units[[CsourceShort]][["biomass"]] = "C-mM"
  
  
  #raw plot
  # viewGroups(dat, yids="biomass", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["biomass"]]= recordPlot()
  
  #smooth biomass
  biomassSM = smoothDat(biomass, dat$Time, nOrder=3, na.rm = T)
  smooths[[CsourceShort]][["biomass"]] = biomassSM$smthObj #save smoothing data
  dat <- addData(dat, "biomassSM", dat=biomassSM$smthDat)
  
  
  #ggplot
  biomassDf = ggplotDf(biomass, concs, times)
  biomassDf = biomassDf[complete.cases(biomassDf),]
  biomassPlt = ggplot(biomassDf)+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y="biomass [C-mM]")
  
  plts[[CsourceShort]][["biomassGG"]] = biomassPlt
  
  
  plts[[CsourceShort]][["items"]][["biomassGG"]] = geom_line(data = biomassDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color="biomass"))
  
  
  biomassDf = data.table(biomassDf)
  biomassDf = biomassDf[,.(value=mean(value)), by=c("groups", "times")]
  biomassDf = data2Perc(biomassDf, settleTime)
  
  plts[[CsourceShort]][["items"]][["biomassStackGG"]] = geom_line(data = biomassDf, aes_(x=quote(times), y=quote(value), group="mean", color="biomass"))
  
  #get yield

  cmol <- dat$layout[["amount"]]
  
  names(cmol) <- dat$layout$well
  
  cmol <- cmol[dat$layout$strain!="blank"] #exclude blanks
  
  SatBeginnings = findSatGrouped(dat, "biomass", Amounts) #find a good beginning for the saturation
  
  if(CSourceRead == "ace"){
    EndTimes = apply(biomass,2,function(x){
      EndTime = which(!is.na(x))
      EndTime = EndTime[length(EndTime)]
    })
    SatTime = SatBeginnings[SatBeginnings<min(EndTimes)]
    SatTime = dat$Time[max(SatTime, na.rm = T)]
  }else{
    SatTime = dat$Time[max(SatBeginnings, na.rm = T)]
  }
  
  bm <- apply(getData(dat, "biomass", xrng=c(SatTime,max(dat$Time))),2,mean, na.rm=T)[names(cmol)] #biomass between highest saturation time and and last time point
  
  
  yieldLm = lm(bm~cmol)
  # if(CsourceShort == "glcSml"){
  #   yieldLm = lm(bm~cmol)
  # }else{
  #   yieldLmFull <- bestFrontLm(x = cmol, y = bm)
  #   yieldLm  = yieldLmFull$lm
  # }
  
  yieldCoef = coef(yieldLm)
  yield = yieldCoef[2]
  
  coefs[[CsourceShort]][["yield"]] = yield #save yield
  
  
  yieldDf = data.frame(cmol, bm)
  
  
  yieldPlt_label = c(sprintf("yield: %.3f", yield), sprintf(": %.3f", summary(yieldLm)$r.squared))
  yieldPlt_label = list(bquote(atop(.(yieldPlt_label[1]),~~~~R^2*.(yieldPlt_label[2]))))
  
  
  yieldPlt = ggplot()+theme_bw()+
    geom_point(data = yieldDf, aes(x=cmol, y=bm))+
    geom_abline(intercept = yieldCoef[1] ,slope=yield, color="red", linetype="dashed")+
    
    labs(x=sprintf("%s concentration [C-mM]",CSource), y="biomass [C-mM]")+
    geom_label(label=yieldPlt_label, aes_(x=max(cmol, na.rm = T),y=0.5*max(bm, na.rm = T)), hjust="right", color="red", parse=T)+
    scale_x_continuous(breaks = concsUniq)+
    expand_limits(y=0)
  
  plts[[CsourceShort]][["yieldGG"]]= yieldPlt
  
  
  #Oxygen
  O2Conc = apply(O2, 2, O2PercToConc) #mmol/l
  
  dat <- addData(dat, "O2Conc", dat=O2Conc)
  units[[CsourceShort]][["O2Conc"]] = "mM"
  
  #save smoothed data
  O2ConcSM = smoothDat(O2Conc, times, na.rm=T, nOrder=3)
  smooths[[CsourceShort]][["O2Conc"]] = O2ConcSM$smthObj
  
  
  #old version of O2 rate
  # # O2AddRate = apply(O2Conc, 2, O2ConcToRate, c_max = O2PercToConc(100), kLa = 260) #mmol/(l h)
  # O2AddRate = (O2PercToConc(100) - O2Conc) * 260 #(c_max - c) * kLa (kLa = Oxygen transfer-coefficient, given by aeration setup) (mmol/(l h)); after Henrys Law
  # 
  # O2Change = (rbind(O2Conc[-1,], NA) - O2Conc) #O2 difference of two consecutive datapoints (mmol/l)
  # row.names(O2Change) = row.names(O2Conc)
  # O2Change = O2Change/ c(dat$Time[-1] - dat$Time[-nTime], NA) #devide by time difference (mmol/(l h))
  # 
  # O2MinusRate = O2Change - O2AddRate # difference in real and theoretical rate is consumption rate (mmol (l h))
  # O2MinusRate = O2MinusRate *10^(-3) # mmol/(ml h) (*1 ml) -> mmol/h
  # 
  # biomassG = biomass/cpg/ mpc * 10^(-6) #g/ml (*1ml) -> g
  # 
  # O2Rate = O2MinusRate/biomassG # O2 consumpion rate of the bacteria (mmol/(h g))
  # 
  # # max(O2Rate[!is.infinite(O2Rate)])
  # # boxplot(O2Rate[!is.infinite(O2Rate)])
  # # 
  # # plot(y = O2Rate, x=matrix(dat$Time, nrow = nrow(O2RateNorm), ncol = ncol(O2RateNorm)))
  # 
  # dat <- addData(dat, "O2Rate", dat=O2Rate)
  # 
  # viewGroups(dat, yids="O2Rate", show.ci95=TRUE, lwd.orig=0)
  # plts[[CsourceShort]][["O2Rate"]]= recordPlot()
  # 
  # 
  # O2RateDf = ggplotDf(O2Rate, concs, times)
  # O2RateDf = O2RateDf[complete.cases(O2RateDf),]
  # 
  # O2RatePlt = ggplot(O2RateDf)+theme_bw()+
  #   geom_line(aes(x=times, y=value, group=as.factor(labels)))+
  #   facet_grid(groups~.)+
  #   geom_hline(yintercept = 0, linetype ="dashed")+
  #   scale_x_continuous(expand=c(0,0))+
  #   labs(x="time [h]", y=expression("norm."~O[2]~"change rate [mmol/(g h)]"))+
  #   bgColsGG+
  #   timeStartGG+
  #   labelFacetGG+
  #   scatterPeaksGG
  # 
  # O2RatePlt  = O2RatePlt + coord_cartesian(ylim = excludedScaling(O2RateDf, excltime = settleTime, exclgroups = 0)) #Only after 0.5h usaful values, so restrict to these values scale, group 0 noisy
  # plts[[CsourceShort]][["O2RateGG"]] = O2RatePlt
  
  
  
  #pH
  pH = getData(dat,"pH")
  
  HPlus = 10^(-pH) #mmol/ml
  
  dat <- addData(dat, "HPlus", dat=HPlus)
  units[[CsourceShort]][["O2Conc"]] = "mmol/ml"
  
  HPlusSM = smoothDat(HPlus, times, na.rm=T, nOrder=3)
  smooths[[CsourceShort]][["HPlus"]] = HPlusSM$smthObj
  
  #old version of HPlus rate
  # HPlusRate = (rbind(HPlus[-1,], NA) - HPlus) #O2 difference of two consecutive datapoints (mmol/ml) (*1ml) -> mmol
  # row.names(HPlusRate) = row.names(HPlus)
  # HPlusRate = HPlusRate/ c(dat$Time[-1] - dat$Time[-nTime], NA) #devide by time difference (mmol/h)
  # HPlusRate = HPlusRate/biomassG #H+ change per g biomass  mmol/(h g)
  # 
  # dat <- addData(dat, "HPlusrate", dat=HPlusRate)
  # 
  # HPlusRateDf = ggplotDf(HPlusRate, concs, times)
  # HPlusRateDf = HPlusRateDf[complete.cases(HPlusRateDf),]
  # 
  # HPlusRatePlt = ggplot(HPlusRateDf)+theme_bw()+
  #   geom_line(aes(x=times, y=value, group=as.factor(labels)))+
  #   facet_grid(groups~.)+
  #   geom_hline(yintercept = 0, linetype ="dashed")+
  #   scale_x_continuous(expand=c(0,0))+
  #   labs(x="time [h]", y=bquote("norm."~H^"+"*~"concentration change rate [mmol/(g h)]"))+
  #   bgColsGG+
  #   timeStartGG+
  #   labelFacetGG+
  #   scatterPeaksGG
  # 
  # HPlusRatePlt = HPlusRatePlt + coord_cartesian(ylim = c(0.001,-0.0005)) #only after 3h useful values
  # plts[[CsourceShort]][["HPlusRateGG"]] = HPlusRatePlt
  
  
  biomassG = biomass/ (1e6 * cpg * mpc) #g/ml
  units[[CsourceShort]][["biomassG"]] = "g/ml"
  
  #get normalized variables and rates
  for(j in list(list("scatter", "norm. scatter alteration rate [AU/(g h)]", "norm. scatter [AU/(g)]"),
                list("ribof", "norm. riboflavin alteration rate [AU/(g h)]", "norm. riboflavin [AU/(g)]"),
                list("O2Conc", bquote("norm."~O[2]~"concentration\nalteration rate [mmol/(l g h)]"),bquote("norm."~O[2]~"[mmol/(l g)]")),
                list("pH","norm. pH alteration rate [1/(g h)]", "norm. pH [1/(g)]"),
                list("NADH", "norm. NADH\n alteration rate [AU/(g h)]", "norm. NADH [AU/(g)]"),
                list("HPlus",bquote("norm."~H^"+"*~"concentration alteration rate [mol/(l g h)]"), bquote("norm."~H^"+"*~"concentration [mol/(l g)]"))
                )
  ){
    varDat = getData(dat, j[[1]])
    varDatNorm = varDat/ biomassG #(X * ml)/g -> X/g
    dat <- addData(dat, paste0(j[[1]],"BMNorm"), dat=varDatNorm)
    
    varSM = smooths[[CsourceShort]][[j[[1]]]]
    
    varSMDeriv = predictSMMat(varSM, times, 1)
    varSMDeriv = varSMDeriv/ biomassG
    
    dat <- addData(dat, paste0(j[[1]],"BMNRate"), dat=varSMDeriv)
    
    #plot normalized rate
    varSMDerivDf = ggplotDf(varSMDeriv, concs, times)
    varSMDerivDf = varSMDerivDf[complete.cases(varSMDerivDf),]
    
    varSMDerivPlt = ggplot(varSMDerivDf)+theme_bw()+
      geom_hline(yintercept = 0, linetype ="dashed")+
      bgColsGG+
      timeStartGG+
      labelFacetGG+
      scatterPeaksGG+
      POIanno+
      
      geom_line(aes(x=times, y=value, group=as.factor(labels)))+
      
      facet_grid(groups~.)+
      scale_x_continuous(expand=c(0,0))+
      labs(x="time [h]", y=j[[2]], color="sample number")
    
    plts[[CsourceShort]][[paste0(j[[1]],"BMNRateGG")]] = varSMDerivPlt
    
    
    #also plot normalized variable
    varDatNormDf = ggplotDf(varDatNorm, concs, times)
    varDatNormDf = varDatNormDf[complete.cases(varDatNormDf),]
    
    varDatNormPlt = ggplot(varDatNormDf)+theme_bw()+
      geom_hline(yintercept = 0, linetype ="dashed")+
      bgColsGG+
      timeStartGG+
      labelFacetGG+
      scatterPeaksGG+
      POIanno+
      
      geom_line(aes(x=times, y=value, group=as.factor(labels)))+
      
      facet_grid(groups~.)+
      scale_x_continuous(expand=c(0,0))+
      labs(x="time [h]", y=j[[3]], color="sample number")
    
    plts[[CsourceShort]][[paste0(j[[1]],"BMNGG")]] = varDatNormPlt
    
    
    #save data as items and exclude blanks
    plts[[CsourceShort]][["items"]][[paste0(j[[1]],"BMNRateGG")]] = geom_line(data = varSMDerivDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color=paste0(j[[1]], "BMNRate")))
    
    plts[[CsourceShort]][["items"]][[paste0(j[[1]],"BMNGG")]] = geom_line(data = varDatNormDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color=paste0(j[[1]], "BMN")))
    
    
    #relative data for comparison
    varSMDerivDf = data.table(varSMDerivDf)
    varSMDerivDf = varSMDerivDf[,.(value=mean(value)), by=c("groups", "times")]
    varSMDerivDf = data2Perc(varSMDerivDf, settleTime)
    
    plts[[CsourceShort]][["items"]][[paste0(j[[1]],"BMNRateStackGG")]] = geom_line(data = varSMDerivDf, aes_(x=quote(times), y=quote(value), group="mean", color=paste0(j[[1]], "BMNRate")))
    
    
    #save relative data
    varDatNormDf = data.table(varDatNormDf)
    varDatNormDf = varDatNormDf[,.(value=mean(value)), by=c("groups", "times")]
    varDatNormDf = data2Perc(varDatNormDf, settleTime)
    
    plts[[CsourceShort]][["items"]][[paste0(j[[1]],"BMNStackGG")]] = geom_line(data = varDatNormDf, aes_(x=quote(times), y=quote(value), group="mean", color=paste0(j[[1]], "BMN")))
  }
  
  for( k in 1:12){
    kNam = c("scatterBMNorm","ribofBMNorm","O2ConcBMNorm","pHBMNorm","NADHBMNorm", "HPlusBMNorm",
      "scatterBMNRate","ribofBMNRate","O2ConcBMNRate","pHBMNRate","NADHBMNRate", "HPlusBMNRate")[k]
    units[[CsourceShort]][[kNam]] = c("AU/g", "AU/g", "mmol/(lg)", "1/g", "AU/g", "mmol/(mlg)",
                                "AU/(gh)", "AU/(gh)", "mmol/(lgh)", "1/(gh)", "AU/(gh)", "mmol/(mlgh)")[k]
  }
    
  
  
  #get O2 consumption for difference of real and theoretical rate of change
  O2AddRate = (O2PercToConc(100) - O2Conc) * kLa #mmol/(lh)
  O2AddRateNorm = O2AddRate/biomassG #mmol/(lhg)
  
  O2ConsumptRate = (getData(dat, "O2ConcBMNRate") - O2AddRateNorm) * 1e-3 #mmol/(lhg) -> mmol/(hg)
  O2ConsumptRate = - O2ConsumptRate #define consumption rate as positive
  
  dat <- addData(dat, "O2ConsumptBMNRate", dat=O2ConsumptRate)
  units[[CsourceShort]][["O2ConsumptBMNRate"]] = "mmol/(hg)"
  
  
  O2ConsumptRateDf = ggplotDf(O2ConsumptRate, concs, times)
  O2ConsumptRateDf = O2ConsumptRateDf[complete.cases(O2ConsumptRateDf),]
  
  plts[[CsourceShort]][["O2ConsumptBMNRateGG"]] = ggplot(O2ConsumptRateDf)+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("norm."~O[2]~"consumption rate"~q[O[2]]~"[mmol/(g h)]"), color="sample number")
  
  
  #save full data as item
  plts[[CsourceShort]][["items"]][["O2ConsumptBMNRateGG"]] = geom_line(data = O2ConsumptRateDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color="O2ConsumptBMNRate"))
  
  
  O2ConsumptRateDf = data.table(O2ConsumptRateDf)
  O2ConsumptRateDf = O2ConsumptRateDf[,.(value=mean(value)), by=c("groups", "times")]
  O2ConsumptRateDf = data2Perc(O2ConsumptRateDf, settleTime)
  
  plts[[CsourceShort]][["items"]][["O2ConsumptBMNRateStackGG"]] = geom_line(data = O2ConsumptRateDf, aes_(x=quote(times), y=quote(value), group="mean", color="O2ConsumptBMNRate"))
  
  
  
  
  #Hplus vs O2
  HPlsO2Df = cbind(plts[[CsourceShort]]$HPlusBMNRateGG$data,
                   O2 = plts[[CsourceShort]]$O2ConsumptBMNRate$data$value)
  
  HPlsO2Plt = ggplot(HPlsO2Df[HPlsO2Df$times >= settleTime,])+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    labelFacetGG+
    
    geom_point(aes(x=O2, y=value))+
    
    facet_grid(groups~.)+
    labs(x=expression("norm."~O[2]~"consumption rate\n"~q[O[2]]~"[mmol/(g h)]"), y=bquote("norm."~H^"+"*~"alteration rate [mmol/(g  h)]"))
  
  plts[[CsourceShort]][["HPlusRate_O2GG"]] = HPlsO2Plt
  
  
  #plot Hplus alteration per O2 alteration
  HPlsO2Df$value = HPlsO2Df$value/HPlsO2Df$O2
  HPlsO2Df = HPlsO2Df[,c("groups", "labels", "times", "value")]
  
  HPlsPerO2Plt = ggplot(HPlsO2Df)+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(y=expression(atop("norm."~H^"+"*~"alteration rate per","norm."~O[2]~"consumption rate"~q[O[2]])), x=bquote("time [h]"))
  
  plts[[CsourceShort]][["HPlusRatePerO2RateGG"]] = HPlsPerO2Plt + coord_cartesian(ylim =excludedScaling(HPlsO2Df, excltime = c(settleTime, 10)))
  
  #save  as items
  plts[[CsourceShort]][["items"]][["HPlusRatePerO2RateGG"]] = geom_line(data = HPlsO2Df, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color="HPlusRatePerO2Rate"))
  
  
  HPlsO2Df = data.table(HPlsO2Df)
  HPlsO2Df = HPlsO2Df[,.(value=mean(value)), by=c("groups", "times")]
  HPlsO2Df = data2Perc(HPlsO2Df, c(settleTime, 10))
  
  plts[[CsourceShort]][["items"]][["HPlusRatePerO2RateStackGG"]] = geom_line(data = HPlsO2Df, aes_(x=quote(times), y=quote(value), group="mean", color="HPlusBMNRate/O2ConsumptBMNRate"))
  
  
  
  #==================================== dpseg analysis ====================================
  
  #growth rates in biomass data
  
  #biomassSM <- apply(biomass, 2, ma, 5)
  
  ## inspect individual well
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # plot(biomass[,id])
  # lines(biomassSM[,id],col=2, type="b", cex=.5)
  
  #dat <- addData(dat, ID="biomassSM", biomassSM)
  
  #can't use dpseg_plate, since NAs are in biomass
  biomassSM = getData(dat, "biomassSM")
  biomassSMNa = is.na(biomassSM)
  
  dpseg = lapply(wells, function(well){
    bmSM = log(biomassSM[,well])
    
    bool = is.na(bmSM)
    bmSM[bool] = bmSM[!bool][sum(!bool)]
    
    dpseg = dpseg::dpseg(x = times, y = bmSM, verb=0, P=.0001)
    
    return(dpseg)
  })
  
  names(dpseg) <- wells
  class(dpseg) <- "dpsegl"
  
  #dpseg <- dpseg_plate(dat, "biomassSM",P=.0001)
  dpsegs[[CsourceShort]][["biomass"]] = dpseg
  
  ## add to plate express object to inspect
  mdat <- addModel(dpseg, dat, ID="biomassSM.dpseg")
  
  ## add growth rates
  mdat <- addModel(dpseg, mdat, ID="biomassSM.mu.dpseg", add.slopes=TRUE, col="#FF0000")
  units[[CsourceShort]][["biomassSM.mu.dpseg"]] = "1/h"
  
  if(any(biomassSMNa)){
    mdat$biomassSM.dpseg$data[biomassSMNa] = NA
    mdat$biomassSM.mu.dpseg$data[biomassSMNa] = NA
  }
  
  # viewPlate(mdat, yids=c("biomassSM","biomassSM.dpseg"), log="y")
  # viewGroups(mdat, yids=c("biomassSM.dpseg"),log="y")
  
  
  ## inspect a single well:
  
  ### NOTE: requires knowledge of the data structures
  ### use RStudio data browser to inspect
  # id <- "A1"
  # if ( nutrient=="ace" ) id <- "E5"
  # 
  # par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
  # plot(mdat$biomassSM$data[,id],log="y",cex=.5)
  # lines(mdat$biomassSM.dpseg$data[,id], lwd=2)
  # par(new=TRUE)
  # plot(mdat$biomass.mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
  # abline(v=dpseg[[id]]$segments[,"start"])
  # abline(v=dpseg[[id]]$segments[,"end"], col=2)
  # axis(4)
  # mtext("growth rate, h-1", 4, par("mgp")[1])
  
  ## inspect to get growth rates
  
  # viewPlate(mdat, yids=c("biomassSM","biomassSM.mu.dpseg"))
  # plts[[CsourceShort]][["growthratesBiomassPlate"]]= recordPlot()
  
  # viewGroups(mdat, yids=c("biomassSM.mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
  #            #ylim=c(-.01,.5),
  #            xlab="time, h", ylab=expression("biomass growth rate"~mu*","~h^-1),
  #            embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
  # abline(h=0, lwd=2)
  # plts[[CsourceShort]][["growthratesBiomassAll"]]= recordPlot()
  
  
  #ggplot
  growthratesBiomassDf = ggplotDf(getData(mdat, "biomassSM.mu.dpseg"), concs, times)
  growthratesBiomassDf = growthratesBiomassDf[complete.cases(growthratesBiomassDf),]
  
  growthratesBiomassPlt = ggplot(growthratesBiomassDf)+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("biomass growth rate "~mu~" ["~h^{-1}~"]"))
  
  plts[[CsourceShort]][["growthratesBiomassGG"]] = growthratesBiomassPlt + coord_cartesian(ylim = excludedScaling(growthratesBiomassDf, excltime = settleTime))
  
  #save growthrates as items
  plts[[CsourceShort]][["items"]][["growthratesBiomassGG"]] = geom_line(data = growthratesBiomassDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color="growthratesBiomass"))
  
  
  growthratesBiomassDf = data.table(growthratesBiomassDf)
  growthratesBiomassDf = growthratesBiomassDf[,.(value=mean(value)), by=c("groups", "times")]
  growthratesBiomassDf = data2Perc(growthratesBiomassDf, settleTime)
  
  plts[[CsourceShort]][["items"]][["growthratesBiomassStackGG"]] = geom_line(data = growthratesBiomassDf, aes_(x=quote(times), y=quote(value), group="mean", color="growthratesBiomass"))
  
  
  
  #growth rates in raw riboflavin data
  
  #ribofSM <- apply(ribof, 2, ma, 5)
  
  ## inspect individual well
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # plot(ribof[,id])
  # lines(ribofSM[,id],col=2, type="b", cex=.5)
  
  #dat <- addData(dat, ID="ribofSM", ribofSM)
  
  
  dpseg <- dpseg_plate(dat, "ribofSM",P=.0001)
  dpsegs[[CsourceShort]][["ribof"]] = dpseg
  
  ## add to plate express object to inspect
  mdat <- addModel(dpseg, mdat, ID="ribofSM.dpseg")
  # viewPlate(mdat, yids=c("ribofSM","ribofSM.dpseg"), log="y")
  # viewGroups(mdat, yids=c("ribofSM.dpseg"),log="y")
  
  ## add growth rates
  mdat <- addModel(dpseg, mdat, ID="ribofSM.mu.dpseg", add.slopes=TRUE, col="#FF0000")
  units[[CsourceShort]][["ribofSM.mu.dpseg"]] = "1/h"
  
  ## inspect a single well:
  
  ### NOTE: requires knowledge of the data structures
  ### use RStudio data browser to inspect
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # 
  # par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
  # plot(mdat$ribofSM$data[,id],log="y",cex=.5)
  # lines(mdat$ribofSM.dpseg$data[,id], lwd=2)
  # par(new=TRUE)
  # plot(mdat$ribofSM.mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
  # abline(v=dpseg[[id]]$segments[,"start"])
  # abline(v=dpseg[[id]]$segments[,"end"], col=2)
  # axis(4)
  # mtext("growth rate, h-1", 4, par("mgp")[1])
  
  ## inspect to get growth rates
  # viewPlate(mdat, yids=c("ribofSM","ribofSM.mu.dpseg"))
  # plts[[CsourceShort]][["growthratesRibofPlate"]]= recordPlot()
  
  
  # viewGroups(mdat, yids=c("ribofSM.mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
  #            #ylim=c(-.01,.1),
  #            xlab="time, h", ylab=expression("riboflavin growth rate"~mu*","~h^-1),
  #            embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
  # abline(h=0, lwd=2)
  # plts[[CsourceShort]][["growthratesRibofAll"]]= recordPlot()
  
  
  #ggplot
  growthratesRibofDf = ggplotDf(getData(mdat, "ribofSM.mu.dpseg"), concs, times)
  growthratesRibofDf = growthratesRibofDf[complete.cases(growthratesRibofDf),]
  
  growthratesRibofPlt = ggplot(growthratesRibofDf)+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("riboflavin growth rate "~mu~" ["~h^{-1}~"]"))
  
  plts[[CsourceShort]][["growthratesRibofGG"]] = growthratesRibofPlt + coord_cartesian(ylim = excludedScaling(growthratesRibofDf, excltime = settleTime))
  
  #save growthrates as items
  plts[[CsourceShort]][["items"]][["growthratesRibofGG"]] = geom_line(data = growthratesRibofDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color="growthratesRibof"))
  
  
  growthratesRibofDf = data.table(growthratesRibofDf)
  growthratesRibofDf = growthratesRibofDf[,.(value=mean(value)), by=c("groups", "times")]
  growthratesRibofDf = data2Perc(growthratesRibofDf, settleTime)
  
  plts[[CsourceShort]][["items"]][["growthratesRibofStackGG"]] = geom_line(data = growthratesRibofDf, aes_(x=quote(times), y=quote(value), group="mean", color="growthratesRibof"))
  
  
  #growth rates in raw scatter data
  
  #scatterSM <- apply(scatter, 2, ma, 5)
  
  ## inspect individual well
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # plot(scatter[,id])
  # lines(scatterSM[,id],col=2, type="b", cex=.5)
  
  #dat <- addData(dat, ID="scatterSM", scatterSM)
  
  dpseg <- dpseg_plate(dat, "scatterSM",P=.0001)
  dpsegs[[CsourceShort]][["scatterSM"]] = dpseg
  
  ## add to plate express object to inspect
  mdat <- addModel(dpseg, mdat, ID="scatterSM.dpseg")
  # viewPlate(mdat, yids=c("scatterSM","scatterSM.dpseg"), log="y")
  # viewGroups(mdat, yids=c("scatterSM.dpseg"),log="y")
  
  ## add growth rates
  mdat <- addModel(dpseg, mdat, ID="scatterSM.mu.dpseg", add.slopes=TRUE, col="#FF0000")
  units[[CsourceShort]][["scatterSM.mu.dpseg"]] = "1/h"
  
  ## inspect a single well:
  
  ### NOTE: requires knowledge of the data structures
  ### use RStudio data browser to inspect
  # id <- "B5"
  # if ( nutrient=="ace" ) id <- "E5"
  # 
  # par(mfcol=c(1,1), mai=c(.5,.5,.1,.5),mgp=c(1.3,.4,0),tcl=-.25)
  # plot(mdat$scatterSM$data[,id],log="y",cex=.5)
  # lines(mdat$scatterSM.dpseg$data[,id], lwd=2)
  # par(new=TRUE)
  # plot(mdat$scatterSM.mu.dpseg$data[,id], axes=FALSE, type="l",ylab=NA)
  # abline(v=dpseg[[id]]$segments[,"start"])
  # abline(v=dpseg[[id]]$segments[,"end"], col=2)
  # axis(4)
  # mtext("growth rate, h-1", 4, par("mgp")[1])
  
  ## inspect to get growth rates
  # viewPlate(mdat, yids=c("scatterSM","scatterSM.mu.dpseg"))
  # plts[[CsourceShort]][["growthratesScatterPlate"]]= recordPlot()
  # 
  # 
  # viewGroups(mdat, yids=c("scatterSM.mu.dpseg"), show.ci95=FALSE, lwd.orig=0,
  #            #ylim=c(-.01,.1),
  #            xlab="time, h", ylab=expression("scatter growth rate"~mu*","~h^-1),
  #            embed=TRUE,no.par=TRUE, g1.legend=FALSE, g2.legend=FALSE)
  # abline(h=0, lwd=2)
  # plts[[CsourceShort]][["growthratesScatterAll"]]= recordPlot()
  
  
  #ggplot
  growthratesScatterDf = ggplotDf(getData(mdat, "scatterSM.mu.dpseg"), concs, times)
  growthratesScatterDf = growthratesScatterDf[complete.cases(growthratesScatterDf),]
  
  growthratesScatterPlt = ggplot(growthratesScatterDf)+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    timeStartGG+
    labelFacetGG+
    scatterPeaksGG+
    POIanno+
    
    geom_line(aes(x=times, y=value, group=as.factor(labels)))+
    
    facet_grid(groups~.)+
    scale_x_continuous(expand=c(0,0))+
    labs(x="time [h]", y=bquote("scatter growth rate "~mu~" ["~h^{-1}~"]"))
  
  plts[[CsourceShort]][["growthratesScatterGG"]] = growthratesScatterPlt + coord_cartesian(ylim = excludedScaling(growthratesScatterDf, excltime = settleTime))
  
  
  #save growthrates as items
  plts[[CsourceShort]][["items"]][["growthratesScatterGG"]] = geom_line(data = growthratesScatterDf, aes_(x=quote(times), y=quote(value), group=quote(as.factor(labels)), color="growthratesScatter"))
  
  
  growthratesScatterDf = data.table(growthratesScatterDf)
  growthratesScatterDf = growthratesScatterDf[,.(value=mean(value)), by=c("groups", "times")]
  growthratesScatterDf = data2Perc(growthratesScatterDf, settleTime)
  
  plts[[CsourceShort]][["items"]][["growthratesScatterStackGG"]] = geom_line(data = growthratesScatterDf, aes_(x=quote(times), y=quote(value), group="mean", color="growthratesScatter"))
  
  ## TODO - platexpress: fix orig colors!
  # colors.fixed <- FALSE
  # if ( colors.fixed ) {
  #   par(mfcol=c(2,1))
  #   viewGroups(mdat, yids=c("scatterSM"), log="y",embed=TRUE,
  #              show.ci95=FALSE, lwd.orig=1)
  #   viewGroups(mdat, yids=c("mu.dpseg"), embed=TRUE, ylim=c(-.01,.22),
  #              g2.legend=FALSE, show.ci95=FALSE, lwd.orig=1)
  #   abline(h=0, lwd=2)
  
  
  O2ConsRate = plts[[CsourceShort]]$items$O2ConsumptBMNRateGG$data
  BMgrowthrates = plts[[CsourceShort]]$items$growthratesBiomassGG$data
  
  
  
  
  growtratesO2Df = join(BMgrowthrates, O2ConsRate, by=c("groups", "labels", "times"))
  names(growtratesO2Df)[c(4,5)] = c("mu", "O2")
  
  growtratesO2Plt = ggplot(growtratesO2Df[growtratesO2Df$times >= settleTime,])+theme_bw()+
    geom_hline(yintercept = 0, linetype ="dashed")+
    #geom_hline(yintercept = 0, linetype ="dashed")+
    bgColsGG+
    labelFacetGG+
    
    geom_point(aes(x=O2, y=mu))+
    
    facet_grid(groups~.)+
    labs(x=expression("norm."~O[2]~"consumption rate"~q[O[2]]~"[mmol/(g h)]"), y=bquote("biomass growth rate "~mu~" ["~h^{-1}~"]"))
  
  plts[[CsourceShort]][["growthrates_O2GG"]] = growtratesO2Plt
  
  #model maximal growthrates against
  #dpsegs[[Amounts]]$biomass
  
  #save data
  dats[[CsourceShort]] <- mdat
  plts[[CsourceShort]]$mods = list("bgColsGG" = bgColsGG,
                                   "timeStartGG" = timeStartGG,
                                   "labelFacetGG" = labelFacetGG,
                                   "scatterPeaksGG" = scatterPeaksGG,
                                   "POIanno" = POIanno)
}


#============================================ save results ============================================ 

#clean up variables
keeps = c("BMass",
          "coefs",
          "dats",
          "dpsegs",
          "plts",
          "smooths",
          "units",
          "biomassEqivalent",
          "settleTime",
          "BMConvs",
          "kLa",
          "INPATH",
          "RESPATH",
          "BMassRibofLm",
          "BMassScatterLm",
          "bestFrontLm",
          "excludedScaling",
          "findDrops",
          "findPlts",
          "findSat",
          "findSatGrouped",
          "ggplotDf",
          "inflateDat",
          "O2PercToConc",
          "predictBM",
          "predictSMMat",
          "savePlts",
          "savePltsPDF",
          "savePltsSingle",
          "smoothDat",
          "stackPlts",
          "data2Perc",
          "label_parsedMixed",
          "updatePackages")

rm(list=ls()[!ls() %in% keeps])
gc()


save.image(file = paste0(INPATH, "/biolector_evalDAT.RData"))

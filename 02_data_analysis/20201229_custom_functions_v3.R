findSat = function(mat, npoints, timeOnCol = T, na.rm=F){
  #' Find the rightbound range of data which seems to saturate
  #' Searches for the rightbound range of data points which minimizes var(points)/n
  #' Uses a heuristic, binary search where the data is repeatedly split into npoint blocks and the best block is evaluated further
  #' 
  #' @param mat timeseries data as a data.frame
  #' @param npoints number of evaluated points per loop in the binary search
  #' @param timeOnCol is the time on the column dimension, i.e. di rows represent the separate samples?
  #' @param na.rm should NA values be removed for evaluation?
  #' 
  #' @return the position of the leftward border, containing the heuristic best saturation range
  
  # Change orientation of data
  # Depending on which dimension contains the time series
  if(!timeOnCol){mat = t(mat)}
  
  endVal = ncol(mat)
  
  # Reshape dataframe
  mat = reshape2::melt(mat)
  mat$Var2 = as.numeric(mat$Var2)
  
  # Note total adjusted variance and stop if not computable
  startVar = var(mat$value, na.rm = na.rm)
  
  if(is.na(startVar)) return(NA)
  
  minVar = startVar/endVal
  minVar_x = 1
  
  # Define total range of points
  currRange = c(1,endVal)
  
  # Define maximal number of steps possible to ensure function stops
  max_steps = (as.integer(endVal/npoints)+1)*2
  for(x in 1:max_steps){
    # Define points evaluated in current binary search repetition
    evalsAll = round(linspace(x1 = currRange[1],
                              x2 = currRange[2],
                              n = npoints + 2)) #evaluated points
    evalsAll = unique(evalsAll)
    evals = evalsAll[c(-1, -length(evalsAll))]
    
    # Compute variace for each evaluated point
    evals_var=
      sapply(evals, function(newStart){
        var(mat[mat$Var2>=newStart,"value"], na.rm = na.rm)
      })
    
    # Compute adjusted variance as var/n
    evals_var = evals_var/(endVal - evals) #reduce accounted variance by number of included points to prefer larger ranges
    
    # Get lowest adjusted variance
    evals_best = which.min(evals_var)
    evals_min = evals_var[evals_best]
    
    # Compare with currently best adjusted variance and note if better, otherwise return current best
    if(evals_min<minVar){
      minVar_x = evals[evals_best]
      minVar = evals_min
      
      currRange = evalsAll[c(evals_best, evals_best+2)]
    }else{
      return(minVar_x)
    }
  }
  return(minVar_x)
}

findSatGrouped = function(dat, ID, group, npoints = 10){
  #' Use the findSat function with a biolector data variable
  #' Assumes the data is grouped by some variable, e.g. used carbon concentration
  #' 
  #' @param dat a biolector data variable
  #' @param ID name of the measurement inside dat to be evaluated
  #' @param group name of the column in the dat layout by which data should be grouped
  #' @param npoints number of evaluated points per loop in the binary search
  #' 
  #' @return A vector with the position of the leftward border, containing the heuristic best saturation range for each unique group level
  
  # Get data from the biolector data variable
  satDat = getData(dat, ID)
  
  # Get grouping variable
  grps = dat$layout[,group]
  
  # Get blank wells
  blnk = dat$layout$strain == "blank"
  
  # Get the unique groups contained in the grouping variable
  grps_uniq = unique(grps)
  grps_uniq = grps_uniq[!is.na(grps_uniq)]
  
  # Apply the findSat function to each groups subset
  # Don't include blanks
  sapply(grps_uniq, function(i){
    
    
    subCols = (grps == i) & !blnk
    subCols = subCols[order(dat$layout$well)]
    subCols[is.na(subCols)] = F
    
    subDat = satDat[,subCols]
    
    findSat(
      subDat,
      npoints,
      timeOnCol = F,
      na.rm = T
    )
  })
}

bestFrontLm = function(y, x){
  #' Find the best linear model containing only leftbound data
  #' 
  #' @param y response variable
  #' @param x independent variable
  #' 
  #' @return a list containing the selected best linear model (lm) and the x values used in its calculation (includedX)

  # Order data in a dataframe
  df = data.frame(x,y)
  
  # Define all possible right borders of x values to evaluate
  groups = sort(unique(x))
  
  # For each x<x_border and corresponding y calculate a linear model
  lms = lapply(groups[-1], function(i){
    lm(y~x, data=df[x<=i,])
  })
  
  # Get the linear model adjusted r squared
  rsqs = sapply(lms, function(i){
    summary(i)$adj.r.squared
  })
  
  # Get the model with the best adjusted r squared
  bestLm = which.max(rsqs)
  
  # Return the best model and the contained x values
  res = list(includedX = groups[1:(bestLm+1)],
             lm = lms[[bestLm]])
  return(res)
}

O2PercToConc = function(perc, p = 21278.25, Temp=310.15){
  #' Convert an O2 percentage in aqueous solution into a molar concentration
  #'
  #' @param perc percent of O2 saturation,
  #' @param p partial pressure of O2 in Pa (21% * 101,325 Pa)
  #' @param Temp temperature in K
  #' 
  #' @return the O2 concentration in mmol l^-1
  
  # Henrys law : c_max = H_cp * p , max gas concentration in mol/m^3
  
  # Define Henry constant at standard conditions and the temperature-adjustment coefficient
  H_cp0 = 1.3 * 10^(-5) # mol/(m^3*Pa)
  dH_cp = 1507.692 # K, mean in Sander et al.
  
  # Get temperature adjusted Henry coefficient
  H_cp = H_cp0 * exp(dH_cp * ( (1/Temp) - (1/298.15) ))
  
  # Get max O2 concentration in mol/l
  c_max = H_cp * p * 10^(-3) 
  
  # Get O2 concentration from max and given percentage
  c = c_max * perc * 10 #*1000 (mol -> mmol), *0.01 (percent -> fraction)
  
  return(c)
}

ggplotDf = function(mat, groups, times){
  #' Turn a matrix or data.frame into a melted data.frame directly usable for ggplot
  #' 
  #' @param mat a matrix or data.frame containing the samples in its columns
  #' @param groups a vector containing the levels of a grouping variable for each sample in mat
  #' @param times a vector containing the time of each row in mat
  #' 
  #' @return a data.frame with columns value: the data points; groups: the respective grouping from groups; labels: the datapointss row name in mat; times: the datapoints time. Can be directly used for ggplot2 plots
  
  # Shape data
  mat = t(mat)
  mat = as.data.frame(mat)
  mat = cbind(groups, labels = row.names(mat), mat)
  
  # Melt data for use with ggplot2
  mat = reshape2::melt(mat, id.vars = c("groups", "labels"), variable.name = "times")
  
  # If lables, times, or groups were turned into factors, turn them back into their original format
  if(is.factor(mat$labels)) mat$labels = levels(mat$labels)[mat$labels]
  if(is.factor(mat$times)) mat$times = times[mat$times]
  
  if(is.factor(mat$groups)){mat$groups = levels(mat$groups)[mat$groups]}
  
  return(mat)
}

savePlts = function(x = "RPlot", type = "pdf", filePath = "", prefix="", width = 10, height = 6){
  #' Saving plots
  #' 
  #' @description Save plots to the current or a specified directory.
  #' Accepts a list of plots or records the currently displayed plot
  #' 
  #' @param x string or list. If x is a string, the current plot will be recorded and saved with prefix_x being the title. If x is a list, all list items must be plots and will be saved with their respective names.
  #' @param type string. Specifies the type of plot, e.g. pdf, png or svg
  #' @param filepath string. Path to the target directory for saving. defaults to the current working directory
  #' @param prefix string. Optional string to be included in the title of the saved plots as prefix_name to tell apart different experiments
  #' @param width integer. Width of the plot. For the same output size compared to pdf and svg plots, the width and height of png plots should be multiplied by 100
  #' @param height integer. Height of the plot. For the same output size compared to pdf and svg plots, the width and height of png plots should be multiplied by 100
  #' 
  #' @examples
  #' savePlts() # saves the currently displayed plot as "Rplot.pdf" to the current working directory
  #' savePlts("scatter", "png" ,"", "HighGA", 1000, 600) # saves the currently displayed plot as "HighGA_scatter.png" to the current working directory
  #' savePlts(plts) # saves all plots in plts as pdfs with their respectibe list names to the to the current working directory 
  
  xType = typeof(x)
  pltFunc = eval(parse(text=type))
  
  if(nchar(prefix)>0){prefix = paste0(prefix,"_")}
  
  
  if(xType == "character"){
    plt = recordPlot()
    
    pltNam = paste0(prefix, x, ".", type)
    
    pltPath = if (nchar(filePath) == 0) {
      pltNam
      
    } else {
      file.path(filePath, pltNam)
    }
    
    pltFunc(pltPath, width = width, height = height)
    print(plt)
    dev.off()
    
  }else if(xType == "list"){
    nams = names(x)
    
    for(nam in nams){
      
      pltNam = paste0(prefix, nam,".", type)
      
      pltPath = if (nchar(filePath) == 0) {
        pltNam
        
      } else {
        file.path(filePath, pltNam)
      }
      
      pltFunc(pltPath, width = width, height = height)
      print(x[[nam]])
      dev.off()
    }
  }else{
    stop("unknown type for x")
  }
}

excludedScaling = function(df, excltime=NULL, exclgroups=c()){
  exclBool = logical(nrow(df))
  
  if(!is.null(excltime)){
    if(length(excltime) ==1){
      exclBool = exclBool | df$times < excltime
    }else if(length(excltime) ==2){
      exclBool = exclBool | df$times < excltime[1] | df$times > excltime[2]
    }
  }
  
  
  if(!length(exclgroups) == 0) exclBool = exclBool | (df$groups %in% exclgroups)
  
  newRange = range(df[!exclBool, "value"])
  
  newRangeexpand = (newRange[2]-newRange[1]) * 0.05
  
  newRange = newRange + c(-1, 1) * newRangeexpand
  
  return(newRange)
}

smoothDat = function(varDat, time, na.rm=F, nOrder=2){
  res = list()
  naVec = rep(NA, nrow(varDat))
  
  varSM = apply(varDat, 2, function(column){
    if(na.rm){
      naBool = !is.na(column)
      
      if(!any(naBool)) return(NA)
      
      column= column[naBool]
      time = time[naBool]
    }
    
    smthDat = sm.spline(x=time, y=column, norder=nOrder)
  })
  
  res$smthObj = varSM #save smoothing object
  
  res$smthDat = sapply(1:length(varSM), function(x){
    naBool = !is.na(varDat[,x])
    
    if(!any(naBool)) return(naVec)
    
    resCol = naVec
    resCol[naBool] = varSM[[x]]$ysmth
    return(resCol)
  })
  
  dimnames(res$smthDat) = dimnames(varDat)
  
  return(res)
}

inflateDat = function(times, fac){
  timediff = times[-1] - times[-length(times)]
  
  newTimeDiff = rep(timediff/fac, each = fac)
  
  newTime = c(times[1], newTimeDiff)
  
  newTime = sapply(2:length(newTime), function(x){
    newTime[x] <<-newTime[x]+newTime[x-1]
  })
  
  newTime[c(rep(F, fac-1), T)] = times[-1]
  
  newTime = c(times[1], newTime)
  
  return(newTime)
}

findDrops = function(mat, times, groups, derivQuantThresh = 0.05, derivQuantThreshExpand = NULL, inclTimeRange = NULL, exclLast = F, combineGroups=F, doPlt = F){
  if(doPlt){
    par(mfrow=c(5,1))
    par(mar=c(1,1,1,1))
  }
  
  matsmooth = smoothDat(mat, times, nOrder = 4)$smthObj
  
  if(exclLast & length(inclTimeRange)<2) times = times[-length(times)]
  
  #Cut off beginning if requested
  if(!is.null(inclTimeRange)){
    if(length(inclTimeRange)==1) inclTimeRange = c(inclTimeRange, Inf)
    
    timeBool = times>=inclTimeRange[1] & times<=inclTimeRange[2]
    times = times[timeBool]
  }
  
  
  
  matDropsList = lapply(matsmooth, function(smth){
    
    smthPrime = predict(smth, times, nderiv=1)
    smth2Prime = predict(smth, times, nderiv=2)
    
    smthPrimeMean = mean(smthPrime)
    smthPrimeSd = sd(smthPrime)
    smth2PrimeMean = mean(smth2Prime)
    smth2PrimeSd = sd(smth2Prime)
    
    derivThresh = qnorm(derivQuantThresh, smthPrimeMean, smthPrimeSd)
    deriv2ThreshDwn = qnorm(derivQuantThresh, smth2PrimeMean, smth2PrimeSd)
    deriv2ThreshUp = qnorm((1-derivQuantThresh), smth2PrimeMean, smth2PrimeSd)
    
    #take points with strong negative slopes
    smthPrimeBool = as.logical(smthPrime<=derivThresh)
    
    if(!any(smthPrimeBool)) return(NULL)
    
    if(is.null(derivQuantThreshExpand)){
      #include points to both sides of these groups, if they have a large second derivative (end points of the drops may have lowered first derivatives)
      smthPrimeBoolSide = (c(F, smthPrimeBool[-length(smthPrimeBool)]) & (smth2Prime >= deriv2ThreshUp)) |
        (c(smthPrimeBool[-1], F) & (smth2Prime <= deriv2ThreshDwn) )
      
      if(doPlt){
        smthPred = predict(smth, times)
        smthPredLow  = smthPredSide = smthPredRest = smthPred
        
        smthPredLow[!smthPrimeBool] = NA
        smthPredSide[!smthPrimeBoolSide] = NA
        smthPredRest[smthPrimeBoolEx|smthPrimeBoolSide] = NA
        
        plot(x=times, y = smthPredRest, ylim=range(smthPred))
        points(x=times, y = smthPredSide, col="blue")
        points(x=times, y = smthPredLow, col="red")
        
        smthPrimeLow = smthPrimeSide = smthPrimeRest = smthPrime
        
        smthPrimeLow[!smthPrimeBool] = NA
        smthPrimeSide[!smthPrimeBoolSide] = NA
        smthPrimeRest[smthPrimeBool|smthPrimeBoolSide] = NA
        
        plot(x=times, y = smthPrimeRest, ylim=range(smthPrime));
        abline(h=derivThresh, col="red")
        points(x=times, y = smthPrimeSide, col="blue")
        points(x=times, y = smthPrimeLow, col="red")
        
        hist(smthPrime, breaks=50)
        abline(v=derivThresh, col="red")
        
        smth2PrimeLow = smth2PrimeSide = smth2PrimeRest = smth2Prime
        
        smth2PrimeLow[!smthPrimeBool] = NA
        smth2PrimeSide[!smthPrimeBoolSide] = NA
        smth2PrimeRest[smthPrimeBool|smthPrimeBoolSide] = NA
        
        plot(x=times, y = smth2PrimeRest, ylim=range(smth2Prime))
        abline(h=deriv2ThreshUp, col="blue")
        abline(h=deriv2ThreshDwn, col="blue")
        points(x=times, y = smth2PrimeSide, col="blue")
        points(x=times, y = smth2PrimeLow, col="red")
        
        hist(smth2Prime, breaks=50)
        abline(v=deriv2ThreshUp, col="blue")
        abline(v=deriv2ThreshDwn, col="blue", main=NULL)
      }
      
      smthPrimeBool = smthPrimeBool | smthPrimeBoolSide
      
      
      #find beginning of consecutive chosen point ranges
      smthGrps = smthPrimeBool + c(1, !(smthPrimeBool[-length(smthPrimeBool)])) * smthPrimeBool
      
      #appoint groups to consecutive chosen points
      grp = 0
      smthGrps = sapply(smthGrps, function(x){
        if(x == 0){
          return(0)
        }else if(x==2){
          grp <<- grp + 1
          return(grp)
        }else if(x==1){
          return(grp)
        }
      })
      
      grps = 1:grp
      
    }else{
      derivThreshEx = qnorm(derivQuantThreshExpand, smthPrimeMean, smthPrimeSd)
      
      #take points with relieved negative slopes of expansion threshold
      smthPrimeBoolEx = as.logical(smthPrime<=derivThreshEx)
      
      smthPrimeBoolSide = (c(F, smthPrimeBoolEx[-length(smthPrimeBoolEx)]) & (smth2Prime >= deriv2ThreshUp)) |
        (c(smthPrimeBoolEx[-1], F) & (smth2Prime <= deriv2ThreshDwn) )
      
      if(doPlt){
        smthPred = predict(smth, times)
        smthPredLow = smthPredEx = smthPredSide = smthPredRest = smthPred
        
        smthPredLow[!smthPrimeBool] = NA
        smthPredEx[!smthPrimeBoolEx] = NA
        smthPredSide[!smthPrimeBoolSide] = NA
        smthPredRest[smthPrimeBoolEx|smthPrimeBoolSide] = NA
        
        plot(x=times, y = smthPredRest, ylim=range(smthPred))
        points(x=times, y = smthPredSide, col="blue")
        points(x=times, y = smthPredEx, col="orange")
        points(x=times, y = smthPredLow, col="red")
        
        smthPrimeLow = smthPrimeEx = smthPrimeSide = smthPrimeRest = smthPrime
        
        smthPrimeLow[!smthPrimeBool] = NA
        smthPrimeEx[!smthPrimeBoolEx] = NA
        smthPrimeSide[!smthPrimeBoolSide] = NA
        smthPrimeRest[smthPrimeBoolEx|smthPrimeBoolSide] = NA
        
        plot(x=times, y = smthPrimeRest, ylim=range(smthPrime));
        abline(h=derivThresh, col="red")
        abline(h=derivThreshEx, col="orange")
        points(x=times, y = smthPrimeSide, col="blue")
        points(x=times, y = smthPrimeEx, col="orange")
        points(x=times, y = smthPrimeLow, col="red")
        
        hist(smthPrime, breaks=50)
        abline(v=derivThresh, col="red")
        abline(v=derivThreshEx, col="orange")
        
        smth2PrimeLow = smth2PrimeEx = smth2PrimeSide = smth2PrimeRest = smth2Prime
        
        smth2PrimeLow[!smthPrimeBool] = NA
        smth2PrimeEx[!smthPrimeBoolEx] = NA
        smth2PrimeSide[!smthPrimeBoolSide] = NA
        smth2PrimeRest[smthPrimeBoolEx|smthPrimeBoolSide] = NA
        
        plot(x=times, y = smth2PrimeRest, ylim=range(smth2Prime))
        abline(h=deriv2ThreshUp, col="blue")
        abline(h=deriv2ThreshDwn, col="blue")
        points(x=times, y = smth2PrimeSide, col="blue")
        points(x=times, y = smth2PrimeEx, col="orange")
        points(x=times, y = smth2PrimeLow, col="red")
        
        hist(smth2Prime, breaks=50)
        abline(v=deriv2ThreshUp, col="blue")
        abline(v=deriv2ThreshDwn, col="blue", main=NULL)
      }
      
      smthPrimeBoolTemp = smthPrimeBool | smthPrimeBoolEx | smthPrimeBoolSide
      
      
      #find beginning of consecutive chosen point ranges
      smthGrps = smthPrimeBoolTemp + c(1, !(smthPrimeBoolTemp[-length(smthPrimeBoolTemp)])) * smthPrimeBoolTemp
      
      #appoint groups to consecutive chosen points
      grp = 0
      smthGrps = sapply(smthGrps, function(x){
        if(x == 0){
          return(0)
        }else if(x==2){
          grp <<- grp + 1
          return(grp)
        }else if(x==1){
          return(grp)
        }
      })
      
      #only use those groups, which include at least one point below the primary threshold
      grps = sort(unique(smthGrps[smthPrimeBool]))
    }
    
    #summarize in matrix
    matDrops = matrix(nrow=length(grps), ncol=4, dimnames = list(NULL,c("t1", "t2", "n1", "n2")))
    
    for(j in 1:length(grps)){
      currGrp = grps[j]
      grpRange = which(smthGrps == currGrp)
      grpRange = grpRange[c(1, length(grpRange))]
      matDrops[j,] = c(times[grpRange], grpRange)
    }
    
    
    if(combineGroups){
      #two groups can be joines, if their combined slope is lower than expected for 10% of datapoints
      combineSlopeThresh = qnorm(0.1, smthPrimeMean, smthPrimeSd)
      combineRsqThresh = 0.5
      
      j=1
      while(j < nrow(matDrops)){
        testRange = matDrops[j,3] : matDrops[j+1,4]
        testDf = data.frame(time=times[testRange], smth=smth$ysmth[testRange,1])
        
        testLm = lm(smth~time, data = testDf)
        
        testLmSlope = coef(testLm)[2]
        testLmRsq = summary(testLm)$r.squared
        
        if(testLmRsq >= combineRsqThresh & testLmSlope <= combineSlopeThresh){
          matDrops[j,c(2,4)] = matDrops[j+1,c(2,4)]
          matDrops = matDrops[-(j+1),]
          next
        }
        
        j=j+1
      }
    }
    
    return(matDrops)
  })
  
  matPeaksDf = data.frame()
  for(z in 1:length(matDropsList)){
    subMat = matDropsList[[z]]
    if(is.null(subMat)) next
    
    subMat = data.frame(subMat)
    
    subMat$groups = groups[z]
    
    subMat$column = z
    
    matPeaksDf = rbind(matPeaksDf, subMat)
  }
  
  if(!is.null(inclTimeRange)){
    matPeaksDf[,c("n1","n2")] = matPeaksDf[,c("n1","n2")] + sum(!timeBool)
  }
  
  if(doPlt) par(mfrow=c(1,1))
  
  return(matPeaksDf)
}

combineDrops = function(found_drops){
  groups = unique(found_drops$groups)
  
  combined_drops = lapply(groups, function(current_group){
    submat = found_drops[
      found_drops$groups==current_group, 
      c("t1", "t2", "n1", "n2", "groups")
    ]
    submat = submat[order(submat$n1), ]
    
    if (nrow(submat)<=1) {
      return(submat)
    } 
    
    peakMat = submat[1, ]
    
    for (sRow in 2:nrow(submat)) {
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
  
  combined_drops = do.call(rbind, combined_drops)
  return(combined_drops)
}

getDropValues <- function(dat, IDs, combined_drops, wells){
  dat_means = lapply(IDs, function(ID){
    rowMeans(getData(dat, ID)[, wells])
  })
  names(dat_means) <- IDs
  
  drop_values = lapply(1:nrow(combined_drops), function(x){
    dropLims = as.numeric(combined_drops[x, c("n1", "n2")])
    
    res = matrix(nrow = length(IDs),
                 ncol = 2,
                 dimnames = list(NULL, c("xmin", "xmax"))
    )
    
    for (i in 1:length(IDs)) {
      res[i,] = dat_means[[i]][dropLims]
    }
    
    res = as.data.frame(res)
    
    res$variable = IDs
    res$drop_no = x
    
    return(res)
  })
  
  drop_values = do.call(rbind, drop_values)
  
  return(drop_values)
}

predictBM = function(dat, model){
  if(is.null(dim(dat))){
    pred = predict(model, newdata = data.frame(value=dat))
    names(pred) = NULL
  } else{
    pred = apply(dat,2, function(x){predict(model, newdata = data.frame(value=x))})
  }
  return(pred)
}

predictSMMat = function(smthObj, times,nderiv=0){
  
  predObj = sapply(smthObj, function(smthdat){
    predict(smthdat, xarg = times, nderiv = nderiv)
  })
  
  return(predObj)
}

stackPlts = function(pltList, category, stacks=NULL, mods=NULL, confinedY = F, xlimit=c(), timeThresh = settleTime){
  if(is.null(stacks)){
    nams = names(pltList[[category]][["items"]])
    nams = nams[grep("[Ss]tack",nams)]
    nams = nams[nams != "stackPlotBasis"]
    nams = sort(nams)
    
    nams2 = names(pltList[[category]][["mods"]])
    nams2 = sort(nams2)
    
    cat("stack items:\n")
    print(nams)
    cat("\nmodifiers:\n")
    print(nams2)
    
    return()
  }
  
  res = pltList[[category]][["items"]][["stackPlotBasis"]]
  
  if(!is.null(mods)){
    if(length(mods) == 1) if(is.na(mods)) mods = names(pltList[[category]][["mods"]])
    
    for(modif in mods){
      res = res + pltList[[category]][["mods"]][[modif]]
    }
  }
  
  maxRange = matrix(nrow=0, ncol=2, dimnames = list(NULL, c("min", "max")))
  ylimit = c()
  
  for(stack in stacks){
    # relCnt = str_count(stack, "/")
    # 
    # if(relCnt == 1){
    #   stack = strsplit(stack, "/")[[1]]
    #   
    #   addPlt <- pltList[[category]][["items"]][[ stack[1] ]]
    #   addPltDiv <- pltList[[category]][["items"]][[ stack[2] ]]
    #   
    #   addPlt$data = merge(addPlt$data, addPltDiv$data, by=c("groups", "times"), sort=F)
    #   addPlt$data = addPlt$data[,.(groups, times, value=value.x/value.y)]
    #   addPlt$data = data2Perc(addPlt$data, timeThresh)
    #   
    #   addPlt$mapping$colour = paste0(addPlt$mapping$colour, "/\n", addPltDiv$mapping$colour)
    #   
    # }else if (relCnt == 0){
    addPlt = pltList[[category]][["items"]][[stack]]
    # }else{
    #   stop("more than on relative '/' not acceptable")
    # }
    
    
    res = res + addPlt
    
    if(confinedY) maxRange = rbind(maxRange, range(addPlt$data$value))
  }
  
  if(confinedY){
    maxRange = data.table(maxRange)
    maxRange = maxRange[,.(min = min(min), max = max(max))]
    maxRange = rbind(maxRange, data.frame(min=-100, max=100))
    
    maxRange = maxRange[,.(min = max(min), max = min(max))]
    maxRangeExpand = (maxRange[,2]-maxRange[,1])*0.05
    maxRange = maxRange[,.(min = min - maxRangeExpand, max = max + maxRangeExpand)]
    
    ylimit = as.numeric(maxRange)
  }
  
  if(confinedY | length(xlimit) == 2){
    res = res + coord_cartesian(ylim = ylimit, xlim = xlimit)
  }
  
  
  return(res)
}

data2Perc = function(mat, timeThresh=0, grouped = F){
  
  if(grouped){
    
    if(length(timeThresh) ==1){
      maxVal = mat[times>=timeThresh, .(max=max(abs(value))), by=groups]
    }else if(length(timeThresh) ==2){
      maxVal = mat[times>=timeThresh[1] & times<=timeThresh[2], .(max=max(abs(value))), by=groups]
    }
    
    res = mat[maxVal, .(groups, times, value = value/i.max* 100) , on="groups"]
    
  }else{
    
    if(length(timeThresh) ==1){
      maxVal = mat[times>=timeThresh, .(max=max(abs(value)))]
    }else if(length(timeThresh) ==2){
      maxVal = mat[times>=timeThresh[1] & times<=timeThresh[2], .(max=max(abs(value)))]
    }
    
    maxVal = as.numeric(maxVal)
    res = mat[, .(groups, times, value = (value/maxVal* 100))]
  }
  
  return(res)
}

subsetGGPlot = function(plotObj, subsetStr){
  dat = data.table(plotObj$data)
  
  dat = dat[eval(parse(text=subsetStr))]
  
  plotObj$data = dat
  return(plotObj)
}

label_parsedMixed = function (labels, multi_line = TRUE) {
  labels <- label_value(labels, multi_line = multi_line)
  if (multi_line) {
    lapply(unname(labels), lapply, function(values) {
      res = try(c(parse(text = as.character(values))), silent = T)
      
      if(class(res) == "try-error") res = as.character(values)
      return(res)
    })
  }
  else {
    lapply(labels, function(values) {
      values <- paste0("list(", values, ")")
      lapply(values, function(expr) c(parse(text = expr)))
    })
  }
}

# Count the number of decimal places
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# Estimate biomass in two growth phases
estimate_two_phases = function(values, times, switchTime, lm1, lm2){
  
  switchBool = times<=switchTime
  switchValue = values[times==switchTime]
  
  coefs1=coef(lm1)
  if(!all(is.na(lm2))){
    coefs2=coef(lm2)
    
    phase1Max = (switchValue * coefs1[2]) + coefs1[1]
    
    res = list((values[switchBool] * coefs1[2]) + coefs1[1],
               ((values[!switchBool] - switchValue)* coefs2[2]) + phase1Max)
    
  }else{
    res = list((values[switchBool] * coefs1[2]) + coefs1[1],
               rep(NA, sum(!switchBool)))
  }
  
  
  res = do.call(c, res)
  return(res)
}

# Create a background colors for a value vector
choose_bg_colors <- function(values, bg_cols){
  
  # First determine, if only one value is present for the factor or if the steps are of equal size
  varied_factor_steps <- 
    values[-1] - values[-length(values)] %>%
    unique()
  
  mono_spaced_bg <- length(values) == 1
  equally_spaced_bg <- length(varied_factor_steps) == 1
  
  # For only one factor level, the one background color can be used
  if (mono_spaced_bg) {
    res = bg_cols
    
  } else if (equally_spaced_bg) {
    
    # For equally spaced levels, use a linear colorRampPalette
    res = colorRampPalette(
      bg_cols
    )(length(values))
    
  }else{
    
    # For unequal spacing adjust the background color number
    # First scale all factor values by multiple of 10 such that there are no decimal places
    factor_amount_Ndec = max(sapply(values, decimalplaces))
    factor_amount_scale = values * 10^factor_amount_Ndec
    
    # Revalue the factor values as positions in a linear scale, starting with 1 for the smalles value
    factor_amount_scale = factor_amount_scale - min(factor_amount_scale) + 1
    
    #Create a linear scaled colorRampPalette and choose the factor level values according to their position
    res = colorRampPalette(bg_cols)(max(factor_amount_scale))
    res = res[factor_amount_scale]
  }
  
  return(res)
}

# Visualize a vector of hex or rgb colors
show_colors <- function(x){
  df <- data.frame(cols = factor(x))
  plt <- ggplot(df, aes(x = cols, y = 1, fill = cols)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = x, guide = F) +
    labs(x = "",y = "")+
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  print(plt)
}

#create a semicricle annotation of points of interest from a dataframe
create_poi_annotation <- function(annotation_df){
  
  if (! is.data.frame(annotation_df)) {
    return(NULL)
  }
  
  annotation_df$group <- as.character(annotation_df$group)
  annotation_df$pos <- as.numeric(paste0(annotation_df$pos,"Inf"))
  annotation_df$vjust <- mapvalues(
    x = annotation_df$pos,
    to = c("bottom","top"),
    from = c(-Inf, Inf)
  )
  
  annotation_pltobj <- list(
    geom_point(
      data = annotation_df,
      aes(x = x, y = pos),
      size = 8,
      shape = 21,
      fill = "grey"
    ),
    geom_text(
      data = annotation_df,
      aes(x = x, y = pos, vjust = vjust, label = lab),
      size = 3.1
    )
  )
  
  return(annotation_pltobj)
}

# Allow for NA in coef function call
coef_na <- function(object, ...){
  if (all(is.na(object)))  {
    return(NA)
  } else{
    coef(object, ...)
  }
}

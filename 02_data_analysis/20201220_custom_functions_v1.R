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

savePlts = function(pltList, filePath = "",width = 1000, height = 600, svgZoom = 1){
  #' Save plots contained in a plts variable for different datasets
  #' Saves in three formats: svg, pdf and png
  #' 
  #' @param pltList a list contining saved plots or plottable ggplot objects in a structure following List$experiment_id$plots or List$experiment_id$subexperiment_id$plots or a mixture of both.
  #' @param filePath path where the plots should be saved
  #' @param width width argument chosen if using the png plotting option. The width used for svg and pdf plotting is chosen as width/100 * svgZoom
  #' @param height height argument chosen if using the png plotting option. The height used for svg and pdf plotting is chosen as height/100 * svgZoom
  #' @param svgZoom factor multiplied to width and height for saved svg and pdf. corresponds to a 1/svgZoom times magnification in these plots.
  
  cats = names(pltList)
  
  for(category in cats){
    nams = names(pltList[[category]])
    nams = nams[!nams %in% c("mods","items")]
    
    mdlBool = grepl("Mdl$",category)
    
    for(nam in nams){
      
      if(!mdlBool){
        pltPath = paste0(filePath, category, "_", nam)
        
        svg(paste0(pltPath, ".svg"), width = width/100 * svgZoom, height = height/100 * svgZoom)
        print(pltList[[category]][[nam]])
        dev.off()
        
        pdf(paste0(pltPath, ".pdf"), width = width/100 * svgZoom, height = height/100 * svgZoom)
        print(pltList[[category]][[nam]])
        dev.off()
        
        png(paste0(pltPath, ".png"), width = width, height = height)
        print(pltList[[category]][[nam]])
        dev.off()
        
      }else{
        subnams = names(pltList[[category]][[nam]])
        
        for(subnam in subnams){
          pltPath = paste0(filePath, category, "_", nam, toupper(sub("^(.).*", "\\1", subnam)), sub("^.(.*)", "\\1", subnam))
          
          svg(paste0(pltPath, ".svg"), width = width/100 * svgZoom, height = height/100 * svgZoom)
          print(pltList[[category]][[nam]][[subnam]])
          dev.off()
          
          pdf(paste0(pltPath, ".pdf"), width = width/100 * svgZoom, height = height/100 * svgZoom)
          print(pltList[[category]][[nam]])
          dev.off()
          
          png(paste0(pltPath, ".png"), width = width, height = height)
          print(pltList[[category]][[nam]][[subnam]])
          dev.off()
        }
      }
    }
  }
}

savePltsSingle = function(pltList, pltVec, filePath = "", width = 1000, height = 600, svgZoom = 1){
  for(pltNam in pltVec){
    
    pltNamSplt = strsplit(pltNam, "_")[[1]]
    category = pltNamSplt[1]
    nam = pltNamSplt[2]
    
    mdlBool = grepl("Mdl$",category)
    pltPath = paste0(filePath, pltNam)
    
    if(!mdlBool){
      svg(paste0(pltPath, ".svg"), width = width/100 * svgZoom, height = height/100 * svgZoom)
      print(pltList[[category]][[nam]])
      dev.off()
      
      pdf(paste0(pltPath, ".pdf"), width = width/100 * svgZoom, height = height/100 * svgZoom)
      print(pltList[[category]][[nam]])
      dev.off()
      
      png(paste0(pltPath, ".png"), width = width, height = height)
      print(pltList[[category]][[nam]])
      dev.off()
      
    }else{
      subnam = regexpr("[A-Z]", nam)
      subnam = substring(nam,subnam)
      subnam = paste0(tolower(sub("^(.).*", "\\1", subnam)), sub("^.(.*)", "\\1", subnam))
      
      nam = substring(nam,1,subnam-1)
      
      svg(paste0(pltPath, ".svg"), width = width/100 * svgZoom, height = height/100 * svgZoom)
      print(pltList[[category]][[nam]][[subnam]])
      dev.off()
      
      pdf(paste0(pltPath, ".pdf"), width = width/100 * svgZoom, height = height/100 * svgZoom)
      print(pltList[[category]][[nam]])
      dev.off()
      
      png(paste0(pltPath, ".png"), width = width, height = height)
      print(pltList[[category]][[nam]][[subnam]])
      dev.off()
    }
  }
}

savePltsPDF = function(pltList, pltVec, pdfName, filePath = "", width = 1000, height = 600){
  pdf(paste0(filePath, pdfName, ".pdf"), width = width, height = height)
  
  for(pltNam in pltVec){
    
    pltNamSplt = strsplit(pltNam, "_")[[1]]
    category = pltNamSplt[1]
    nam = pltNamSplt[2]
    
    mdlBool = grepl("Mdl$",category)
    
    if(!mdlBool){
      print(pltList[[category]][[nam]])
    }else{
      subnamPos = regexpr("[A-Z]", nam)
      subnam = substring(nam,subnamPos)
      subnam = paste0(tolower(sub("^(.).*", "\\1", subnam)), sub("^.(.*)", "\\1", subnam))
      
      nam = substring(nam,1,subnamPos-1)
      
      print(pltList[[category]][[nam]][[subnam]])
    }
  }
  dev.off()
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
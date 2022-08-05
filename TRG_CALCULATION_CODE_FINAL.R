library(runner)
library(ggplot2)
library(tidyverse)
library(car)
library(ggrepel)
library(gridExtra)
library(factoextra)
library(cowplot)
library(stringr)

setwd("~/Documents/wale_docs/phd/data") # optional

#df <- read.table(file = "test_ground/Output/Ed1a_10.txt", sep = '\t', header = T)

#Auxilliary Functions
#!important to run the derivative aux function first
derivative <- function(input, pl=F, adjParam = 0.1){#you can set pl = False to reduce the number of plots generated in the main function
  
  minData <- min(input$wellName, na.rm = T)
  #print(minData)
  
  input <- subset(input, wellName < minData+adjParam)
  #print(length(input$Time))
  #print(input$wellName)
  grads= c()
  slope <- runner(
    x = input,
    k = 4,
    idx =  c(1:length(input$Time)),
    function(x) {
      coefficients(lm(wellName ~ Time, data = x))[[2]]
    }
  )
  
  grads <- c(grads, slope)
  #print(grads)
  
  
  roundGrads <- round(grads, digit=3)
  #print(roundGrads)
  maxGrads = NULL
  
  if (length(roundGrads) > 3) {
    selection <- roundGrads[4:length(roundGrads)]
    maxGrads = which(roundGrads == max(selection, na.rm = T))-4
  }else{
    maxGrads = which(roundGrads == max(roundGrads, na.rm = T))-4
  }
  
  
  
  secondInput <- data.frame('Time' = input$Time, 'grads'=grads)
  #print()
  
  secondGrads= c()
  slope2 <- runner(
    x = na.omit(secondInput),
    k = 4,
    #idx =  c(1:length(secondInput$grads)),
    function(x) {
      coefficients(lm(grads ~ Time, data = x))[[2]]
    }
  )
  
  secondGrads <- c(secondGrads, slope2)
  #print(secondGrads)
  #6
  #print(tail(secondGrads, 35))
  
  if (length(subset(input, wellName > minData, wellName < minData+0.2)$Time) > 6) {
    solution <- which(secondGrads == max(tail(secondGrads[5:(length(secondGrads)-10)], -10), na.rm = T))
  } else {
    solution <- which(secondGrads == max(tail(secondGrads, 6), na.rm = T))
  }
  
  
  OD <- NULL
  
  #print(solution)
  if (is.finite(solution) && length(solution) == 1) {
    OD <- input$wellName[solution]
    if (pl==T) {
      #plot(input$wellName, ylim = c(-0.2,0.4))
      #lines(secondGrads, col='blue')
      #abline(v = solution)
      #abline(h= OD)
    }
  } else {
    solution = 0
    OD = 0
  }
  
  
  return(c(solution, OD,maxGrads))
}


#NOTE
#TRG_CALC can function on its own if you supply an already imported dataset i.e. the read.table() importer

TRG_CALCULATOR_THREE <- function(plate, windowSize=4, date, plate_no, method, colorCode='chocolate', adjP = 0.1, od_transf_start = -4) {
  
  
  xAxis = plate$Time
  wellName = c(colnames(plate[2:length(plate)]))
  plateTRG = NULL
  growthRate = NULL
  slopeDetails = NULL
  corData = NULL
  RSQ = NULL
  pbl = NULL
  distance = NULL
  elbowTime = NULL
  winDiff = NULL
  INTERCEPT = NULL
  
  for (column in 2:ncol(plate)){
    print(colnames(plate[column]))
    cd = plate[,column]
    cdDist = plate[,column]
    
    #determine baseline
    baseline_dataframe <- data.frame('wellName' = plate[,column][3:nrow(plate)], 'Time' = plate[,'Time'][3:nrow(plate)])
    baseline_range <- derivative(baseline_dataframe, adjParam = adjP)[1]
    #print(baseline_range)
    if (baseline_range>7) {
      ranger <- c(baseline_range-1, baseline_range-2, baseline_range-3)
    }else{
      ranger <- c(baseline_range, baseline_range+1)
    }
    #print(ranger)
    baseline <- mean(cd[ranger[ranger>0]], na.rm = T)
    elbowPoint <- head(baseline_dataframe$Time[baseline_range-3], n = 1) #track forward to select true position in the unedited data
    baseElbowPoint <- cd[baseline_range]
    #print(baseline)
    #print(paste('elbowPoint =', elbowPoint))
    searchBase = cd[baseline_range]
    
    distSol <- 0
    
    #remove no growth baselines
    if (baseline_range>0) {
      baseline_range <- baseline_range
    }else {
      baseline_range <- length(plate[,'Time'])
    }
    
    if (baseline_range>0) {
      baseline <- baseline
    }else {
      baseline <- cd[length(plate[,'Time'])]
    }
    
    
    if (!is.na(tail(elbowPoint[1]))) {
      elbowPoint <- elbowPoint
    }else {
      elbowPoint <- xAxis[length(plate[,'Time'])]
    }
    
    plot(xAxis, cd, pch=1, col=colorCode, ylim = c(0,1.2), xlim=c(0,24), 
         main = colnames(plate[column]), cex=0.5, ylab = 'Raw', xlab = 'Time')
    abline(v = elbowPoint, col = 'grey')
    abline(h = baseline, col= 'blue')
    
    #print(baseline_range)
    
    pickup <- head(baseline_range, n = 1)
    if (pickup < 1) {
      pickup = 1
    }else {
      pickup = pickup
    }
    cd = plate[,column][pickup:nrow(plate)]
    xData = plate[,'Time'][pickup:nrow(plate)]
    yAxis <- cd-baseline
    utcorData <- data.frame('x' = xData, 'y'= suppressWarnings(log(yAxis)))
    
    corData <- subset(utcorData, y>=log(exp(od_transf_start)) & y<=log(0.5))
    
    if (nrow(corData) > 1) {
      #run sliding window
      slope <- runner(
        x = corData,
        k = windowSize,
        at = c(3:windowSize),
        function(x) {
          coefficients(lm(y ~ x, data = x))[2]
        }
      )
      # save index of the max slope
      whichSlope <- c(which(slope == max(slope, na.rm = T)))[[1]] 
      # save the max Slope == growth rate
      trueMaxSlope <- slope[whichSlope][[1]]
      #print(trueMaxSlope)
      
      slopeUsed <- runner(
        x = corData,
        k = 3,
        at = c(3,4),
        function(x) {
          lm(y ~ x, data = x)
        }
      )
      
      #print('wale')
      #print(slopeUsed)
      #print(whichSlope)
      
      slopeDetails <- slopeUsed[,whichSlope]
      transRawD <- (slopeUsed[,whichSlope]$model)$y
      indTransRawD = which(utcorData$y %in% transRawD)+ pickup
      fitness <- cor(x = slopeUsed[,whichSlope]$model$y, 
                     y = slopeUsed[,whichSlope]$fitted.values, 
                     method = 'pearson', use = "complete.obs")
      #min x used in calculation
      rateMin <- min((slopeUsed[,whichSlope])$model$y)
      #max y used in calculation
      rateMax <- max((slopeUsed[,whichSlope])$model$y)
      #intercept
      intercept <- coefficients(slopeDetails)[[1]]
      #print(paste('intercept = ', intercept))
      #slope 
      slope <- coefficients(slopeDetails)[[2]]
      #print(slope)
      #TRG Calculator USING new n= 0.01123556
      TRG_absolute <- (log(0.0035) - intercept)/slope
      
      #print(TRG_absolute)
      
      if (TRG_absolute < 0 | slope < 0) {
        TRG_absolute <- NA
        slope <- NA
      }else{
        TRG_absolute <- TRG_absolute
        slope <- slope
      }
      
      #plot point used for fitting in raw data
      points(x = plate[,'Time'][indTransRawD] , y = plate[,column][indTransRawD], 
             col = 'red', pch = 19)
      
      wDiff <- min(plate[,column][indTransRawD], na.rm =T) - baseline
      
      #BASELINE DISTANCE CALCULATION
      pBaseline <- rep(baseline, baseline_range)
      qPoints <- cdDist[c(1:baseline_range)]
      #print(pBaseline)
      #print(qPoints)
      xDist <- rbind(pBaseline,qPoints)
      distSol = stats::dist(xDist, method = "manhattan")
      
      #plot fits
      plot(corData, ylim = c(-6,0), xlim = c(0,24), main =colnames(plate[column]), cex=0.75 )
      points(x = (slopeUsed[,whichSlope]$model)$x, y = (slopeUsed[,whichSlope]$model)$y, col = 'red', pch=19)
      abline(slopeDetails)
      abline(h = log(0.03019738))
      abline(h = log(0.4))
      text(x = 0, y = 0, labels = paste('GR: ', round(slope, 2)), cex =1, pos = 4)
      text(x = 0, y = -0.5, labels = paste('TRG: ', round(TRG_absolute, 2)), cex =1, pos = 4)
      text(x = 0, y = -1.0, labels = paste('RSQ: ', round(fitness, 7)), cex =1, pos = 4)
      text(x = 0, y = -1.5, labels = paste('dist: ', c(distSol)), 7, cex =1, pos = 4)
    }else{
      TRG_absolute <- NA
      slope <- NA
      Rsq_fit <- 0
      fitness <- 0
      plot(0,24, ylim = c(-6,0), xlim = c(0,24), main =colnames(plate[column]))
      text(x = 7.5, y = -3.0, labels = paste('No growth'), cex =1, pos = 4)
      fitness <- NA
      distSol <- NA
      elbowPoint <- NA
      wDiff <- NA
      intercept <- NA
    }
    
    plateTRG <- rbind(plateTRG, TRG_absolute)
    growthRate <- rbind(growthRate, slope)
    ID <- paste(date, plate_no, wellName, sep = '_')
    RSQ <- rbind(RSQ, fitness)
    pbl <- rbind(pbl, baseline)
    distance <- rbind(distance, distSol)
    elbowTime <- rbind(elbowTime, elbowPoint)
    winDiff <- rbind(winDiff, wDiff)
    INTERCEPT <- rbind(INTERCEPT, intercept)
  }
  
  plateReport <- data.frame('wellName' = wellName, 
                            'TRG' = plateTRG, 
                            'GRate'= growthRate,
                            'Identifier' = ID,
                            'RSQ' = RSQ,
                            'pbl' = pbl,
                            'distance' = distance,
                            'elbowTime' = elbowTime, 
                            'winDiff' = winDiff,
                            'Y_INTERCEPT' = INTERCEPT)
  
  return(plateReport)
}

#a = TRG_CALCULATOR_THREE(plate = litData, date = NULL, plate_no = NULL, min_start = -7)


#NOTE
#TO calculate for multiple datasets in multiple folders, use plottingMachine. 
#IT WORKS WITH THE LAB'S CODING FOLDER STRUcture, "../Output" 

#Ploting function
#before using it, you can define create a pdf file to dump the plots !optional
plottingMachineX1 <- function(folder_Name){
  pwd <- '~/Documents/wale_docs/phd/data/' # a superfolder with folders containing a subfolder named 'Output' !important
  
  All_files = list.files(paste0(pwd, folder_Name, '/Output/', collapse = NULL))
  Plates = All_files[grep("txt",All_files)]
  
  #If error, use print(Plates) to debug
  #print(Plates)
  
  par(mfrow = c(1,1))
  plot(0,0, col='white', axes = F)
  text(x = 0, y = 0, folder_Name, cex=3)
  comprehensive_solution <- NULL
  
  for (plate in Plates) {
    newPWD <- paste0(pwd, folder_Name, '/Output/',plate)
    print(newPWD)
    plate_data <- read.table(file = newPWD, sep = '\t', header = T)
    print(paste('Plotting ', plate))
    par(mfrow = c(1,1))
    plot(0,0, col='white', axes = F)
    text(x = 0, y = 0, plate, cex=3)
    par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
    
    identifier_code = str_replace_all(folder_Name, '_', '')
    
    solution = NULL
    if (length(plate_data) == 97) {
      solution = TRG_CALCULATOR_THREE(plate_data, date = identifier_code, plate_no = which(plate == Plates))
    }
    
    comprehensive_solution <- rbind(comprehensive_solution, solution)
  }
  
  return(comprehensive_solution) 
}

#Direct Implementation
sample_data <- read.table(file = "/Users/josoga2/Documents/wale_docs/phd/data/RG_JV_20h_96wp_220729_ref.txt", sep = '\t', header = T)
#litData <- data.frame('Time' = p1$Time, 'Well' = p1$O23)
pdf(file = 'mystery_data_FIXED.pdf' ,useDingbats = F, paper = 'a4r', width = 120, height = 80)
par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
a = TRG_CALCULATOR_THREE(plate = sample_data, date = NULL, plate_no = NULL)

dev.off()


#implementation: Just specify a folder with output in it
#you can also specify multiple folders with the vector list structure c(folder_1, folder_2, ..., folder_n)

#I recommend this pdf size, for reasons I can't explain. It just works
#pdf(file = 'fits_ALL.pdf',useDingbats = F, paper = 'a4r', width = 120, height = 80) 

#result <- plottingMachineX1(folder_Name = 'test_ground')
#export your data
#write.table(x = result, file = 'simple_Result.csv', quote = T, sep = ',',row.names = F)

#dev.off()

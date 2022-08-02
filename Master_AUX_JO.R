#source
setwd("/Users/josoga2/Documents/wale_docs/phd/data")
source("/Users/josoga2/Documents/wale_docs/phd/R/Collection_of_scripts/Aux_functions.R")

library(runner)
library(ggplot2)
library(tidyverse)
library(car)
library(ggrepel)
library(gridExtra)
library(factoextra)
library(cowplot)
library(stringr)
library(bioassays)
library(ggplot2)
library(phenoScreen)
library(gridExtra)
library(ggpubr)
library(car)
library(ggvenn)
library(RColorBrewer)
library(reshape2)


#quadrant_splitter

wellPlatePlotter <- function(SingleColData, plateType, labelTitle) {
  
  if (plateType == 96) {
    mymeta <- metafile96
    mymeta$type <- 'S1'
    mymeta$id <- 'Sample'
    mymeta$concentration <- 0.15
    rawdat <- cbind('X' = rawdata96$X , t(data.frame(matrix(data = SingleColData, nrow = 12, ncol = 8) )))
    #print(rawdat)
    rawdat <- data2plateformat(data = rawdat, platetype = 96)
    trg_df<- plate2df(rawdat)
    
    pll <- ggplot(data = trg_df, aes(x = reorder(x = row, desc(row)), y =col, col=value))+
      geom_point(size=18)+ylim(c(0,13))+coord_flip()+theme_bw() +
      scale_y_continuous(position = 'right' ,breaks=c(1:12), labels=c(1:12))+ 
      scale_color_gradient2(low = '#2a6a99', mid = '#c6dbef', midpoint = 10, high = '#d88546', limits= c(0,20))+
      labs(x=NULL, y = NULL)+ theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_text(label= round(c(trg_df$value), 2), color = 'black', size=5)+
      ggtitle(label = paste0('Plate '))
  }
  
  if (plateType == 384) {
    mymeta <- metafile384
    mymeta$type <- 'S1'
    mymeta$id <- 'Sample'
    mymeta$concentration <- 0.15
    rawdat <- cbind('X' = rawdata384$X , t(data.frame(matrix(data = SingleColData, nrow = 24, ncol = 16) )))
    #print(rawdat)
    rawdat <- data2plateformat(data = rawdat, platetype = 384)
    trg_df<- plate2df(rawdat)
    
    pll <- ggplot(data = trg_df, aes(x = reorder(x = row, desc(row)), y =col, col=value))+
      geom_point(size=12.5)+ylim(c(0,13))+coord_flip()+theme_bw() +
      scale_y_continuous(position = 'right' ,breaks=c(1:24), labels=c(1:24))+ 
      scale_color_gradient2(low = '#2a6a99', mid = '#c6dbef', midpoint = 10, high = '#d88546', limits= c(0,20))+
      labs(x=NULL, y = NULL)+ theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_text(label= round(c(trg_df$value), 2), color = 'black', size=3)+
      ggtitle(label = paste0( labelTitle))
  }
  return(pll)
}


serialPlotter <- function(whereIsFileLocated, 
                          saveFileAs, 
                          plateType,
                          sepStat = 'tog',
                          Treatment = 'Treatment') { 
  
  All_files = list.files(paste0(whereIsFileLocated, collapse = NULL))
  Plates = All_files[grep("txt",All_files)]
  n <- 0 
  my_curr_plates <- c()
  
  for (plate in Plates) {
    n <- n+1
    my_curr_ <- read.table(paste0(whereIsFileLocated, plate), header = T)
    my_curr_plates <- c(my_curr_plates, list(my_curr_))
  }
  
  
  pdf(file = paste0(whereIsFileLocated, saveFileAs) ,useDingbats = F, paper = 'a4r', width = 120, height = 80)
  #Replicate Correlation
  par(mfrow = c(1,1))
  plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
  text(0,0, 'Replicate Correlation')
  
  n <- 1
  colIDs <- c()
  melted_mtx <- data.frame('Plate1' = melt(my_curr_plates[[1]][2:length(my_curr_plates[[1]])])$value)
  for (currNo in 2:length(my_curr_plates)) {
    n <- n+1
    colIDs  <- c(colIDs, paste0('Plate', n))
    curr <- my_curr_plates[[currNo]]
    newCol <- melt(curr[2:length(curr)])$value
    melted_mtx[colIDs] <- newCol
  }
  
  pairs(melted_mtx, col = 'black', cex= 0.75, pch=21, bg='lightgreen', upper.panel = panel.cor, lower.panel = panel.smooth)
  
  ####create raw plots
  par(mfrow = c(1,1))
  plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
  text(0,0, 'Experiment 1')
  
  par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
  wellNames <- names(my_curr_plates[[1]])[2:length(my_curr_plates[[1]])]
  
  if (sepStat == 'tog') {
    for (well in wellNames) {
      col_n <- 0
      plot(0,0, xlim = c(0,20), ylim = c(0,1.2), col='white', main = well)
      legend('topleft', legend = c('unt', 'rep1','rep2','rep3'), col = 1:4, pch = 19, cex = 0.75)
      for (plate in my_curr_plates) {
        col_n <- col_n+1
        points(plate[['Time']], plate[[well]], col = col_n)
      }
    }
  } else {
    for (plate in my_curr_plates) {
      col_n <- 0
      for (well in wellNames) {
        plot(plate[['Time']], plate[[well]], col = 4, pch = 19, cex = 0.75, 
             main = well, xlim = c(0,20), ylim = c(0,1.2))
      }
    }
  }
  
  
  
  # compute time of regrowth for each plate completely
  all_Data_HOLDER <- c()
  
  par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
  for (plate in my_curr_plates) {
    TEMP_DATA_HOLDER <- TRG_CALCULATOR_THREE(plate = plate, date = NULL, 
                                                 plate_no = paste0('Plate ', n), 
                                                 adjP = 0.15)
    
    all_Data_HOLDER <- c(all_Data_HOLDER, list(TEMP_DATA_HOLDER))
  }
  
  #take only TRG for plotting
  TRG_HOLDER <- c()
  for (holding in all_Data_HOLDER) {
    TRG_HOLDER <- c(TRG_HOLDER, list(holding$TRG))
  }
  
  #take only GRate for plotting
  GRATE_HOLDER <- c()
  for (holding in all_Data_HOLDER) {
    GRATE_HOLDER <- c(GRATE_HOLDER, list(holding$GRate))
  }
  
  #PLOT PLATES WITH PLATEPLOTTER
  N <- 0
  for (compTRG in TRG_HOLDER) {
    N = N+1
    ppll <- wellPlatePlotter(SingleColData = c(compTRG), 
                            plateType = plateType, labelTitle = paste0('plate', N))
    
    grid.arrange(ppll)
  }
  print('COMPLETED WELLPLATES')
  
  #plot TRG VS GRATE
  NPLOTS <- list()
  N <- 0
  for (comp in all_Data_HOLDER) {
    N <- N+1
    NPLOTS[[N]] <- ggplot(comp, aes(x = TRG, y = GRate, color = RSQ))+
      geom_point()+ggtitle(paste0('Plate ', N))+
      scale_color_gradient(low = '#24749B', high = '#F6992F')+
      lims(x = c(0,20), y = c(0,2.5))
  }
  
  do.call(grid.arrange,NPLOTS)
  
  
  #outliers with replicates
  par(mfrow = c(round(length(TRG_HOLDER)/2), 2))
  OUTLIER_HOLDER <- c()
  for (holding in TRG_HOLDER) {
    OUTLIER_HOLDER <- c(OUTLIER_HOLDER, list(Boxplot(holding, notch =T, col =8)))
  }
  
  OUTLIER_LABEL_HOLDER <- c()
  for (plateOrder in 1:length(all_Data_HOLDER)) {
    OUTLIER_LABEL_HOLDER <- c(OUTLIER_LABEL_HOLDER, 
                              list(all_Data_HOLDER[[plateOrder]][OUTLIER_HOLDER[[plateOrder]],]$wellName))
  }
  
  #PLOT VENN DIagram for intersection of outliers
  if (sepStat == 'tog') {
    veggven <- ggvenn(data = OUTLIER_LABEL_HOLDER[2:4], columns = 1:3, stroke_size = 0,
                      fill_color = c('#F6992F', '#24749B', '#22B258'), fill_alpha = 0.5, text_size = 5, show_elements = F)
    grid.arrange(veggven)
  }
  
  #boxplot to emphasize outliers
  print('completed ggven')
  
  if (sepStat == 'tog') {
    par(mfrow = c(1,1), mar = c(5, 5, 5, 5))
    
    #boxplot of all data
    bxpdat <- boxplot(TRG_HOLDER, ylim= c(0,20), outline = F, varwidth = T, col = 'grey', names = paste0('Plate ', 1:length(TRG_HOLDER)),
                      boxwex=.8, notch = T, border = 1, lwd=2, cex.lab = 1.5, cex.axis = 1.2,
                      ylab = substitute(paste(bold('Time of Regrowth (Hrs)'))), 
                      xlab=substitute(paste(bold(Treatment))))
    print('completed box plots')
    
    n <- 0
    for (labelling in 1:length(OUTLIER_LABEL_HOLDER)) {
      bxpdat <- boxplot(TRG_HOLDER[[1]], plot = F)
      n <- n+1
      text(x = labelling, y= bxpdat$out, labels = c(OUTLIER_LABEL_HOLDER[[labelling]], ''), col = 1:10, cex=1.5) 
    }
  }
  
  dev.off()
  return(all_Data_HOLDER)
  
}




setwd("~/Documents/wale_docs/phd/data") # optional

#df <- read.table(file = "test_ground/Output/Ed1a_10.txt", sep = '\t', header = T)

#Auxilliary Functions
#!important to run the derivative aux function first
derivative <- function(input, pl=T, adjParam = 0.1){#you can set pl = False to reduce the number of plots generated in the main function
  
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

TRG_CALCULATOR_THREE <- function(plate, windowSize=4, date, plate_no, method, colorCode='chocolate', adjP = 0.1, min_start = -4) {
  
  
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
    
    print(baseline_range)
    
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
    
    corData <- subset(utcorData, y>=log(exp(min_start)) & y<=log(0.5))
    
    if (nrow(corData) > 1) {
      #run sliding window
      slope <- runner(
        x = corData,
        k = 3,
        at = c(3,4),
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
      TRG_absolute <- (log(0.01113222) - intercept)/slope
      
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


#import and analyze plates from SAMI
require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(LSD)
require(corrplot)
library(igraph)
require(plotrix)
require(gtools)
library(bioassays)

#======== Setup path and directories ============
setwd("/Users/josoga2/Documents/wale_docs/phd/")
here_path = "/Users/josoga2/Documents/wale_docs/phd/"
source(paste0(here_path,"/R/Collection_of_scripts/Aux_functions.R",collapse=NULL))

Load_dir = paste0(here_path,"data/Phase_1_Data/26_07_22/") 
Out_dir = paste0(here_path,"data/Phase_1_Data/26_07_22/Output/")

All_files = list.files(Load_dir, pattern = 'log')
#All_files = All_files[[1]] #de-comment this in case you don't want to convert all files int he folder

#======== The conversion function & conversion ============
#Authored by ARB
#MODIFIED BY ME

Ext_SAMIplates <- function(file_path, output_dir)
{
  dataset <- read.table(file_path, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, skip = 4)
  
  #Remove NA columns in case of 96 weel plates recorded as 384
  if(length(dataset)>120 && sum(is.na(dataset[length(dataset)])) == length(dataset[[1]]))
  {
    dataset = dataset[1:(3+96)]
  }
  
  header =c("Family","Plate","Time",paste0(rep(LETTERS[1:8],each=12),seq(1,12)))
  if(length(dataset)>120)
  {header=c("Family","Plate","Time",paste0(rep(LETTERS[1:16],each=24),seq(1,24)))}
  
  names(dataset) = header        
  
  family_vec <- as.numeric(as.vector(dataset$Family))
  nr_families <- max(family_vec)
  
  families_OD_list <- list()
  for (family in 1:nr_families)
  {
    family_OD_data = dataset[grep_exact(family,family_vec),]
    
    if(length(grep(T,duplicated(family_OD_data))>0))
    {family_OD_data = family_OD_data[grep(F,duplicated(family_OD_data)),]}
    
    families_OD_list[[length(families_OD_list)+1]] <- family_OD_data
  }
  
  #Covert Date and time columns to time-intervals in hours and write the files
  nx <- 1
  
  for(family in 1:nr_families)
  {
    family_data = families_OD_list[[family]]
    Date_time <- as.character(as.vector(family_data$Time))
    
    dtm <- strptime(Date_time, format = "%m/%d/%Y %H:%M:%S", tz = "CET")
    Time_h = vector(mode = "numeric", length = length(dtm))  
    for (i in 2:length(Time_h))
    {
      j = i-1
      diff = dtm[i]-dtm[j]
      Time_h[i] = Time_h[j] + as.numeric(diff, units="hours")
    }
    
    family_data$Time = Time_h
    PlateID = as.character(family_data$Plate[1])
    #print(family_data$Plate[1])
    
    family_data = family_data[-c(1,2)]
    
    write.table(family_data,file = paste0(Out_dir,'Plate', nx,  ".txt"), append = F, sep = "\t",row.names = F,quote = F)
    nx <- nx+1
  }
}



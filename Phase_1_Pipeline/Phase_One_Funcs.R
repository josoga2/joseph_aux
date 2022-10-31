#all functions so farS
#source
setwd("/Users/josoga2/Documents/wale_docs/phd/data")
source("/Users/josoga2/Documents/wale_docs/phd/R/Collection_of_scripts/Aux_functions.R")

LibMap <- read.table('../PLANS/LibMap.txt', sep = '\t', quote = '\t', header = T)
LibMap$compoundID <- paste0('LP',LibMap$Library.plate,'_', LibMap$Well)
cleanLibMap <- LibMap[c('compoundID', 'Catalog.Number')]


######## DMSO4 DATA #####
#pick all DMSO wells
allDmsoWell <- unique(subset(LibMap, Catalog.Number == 'DMSO')$Well)
#########################

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

wellPlatePlotter <- function(SingleColData, plateType, labelTitle, colorScale = c(0,25)) {
  
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
      scale_color_gradient2(low = '#2a6a99', mid = '#c6dbef', midpoint = mean(colorScale), high = '#d88546', limits= colorScale)+
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
      scale_color_gradient2(low = '#2a6a99', mid = '#c6dbef', midpoint = mean(colorScale), high = '#d88546', limits= colorScale)+
      labs(x=NULL, y = NULL)+ theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_text(label= round(c(trg_df$value), 2), color = 'black', size=3)+
      ggtitle(label = paste0( labelTitle))
  }
  return(pll)
}


zScorer <- function(plateDataset) {
  std_rep <- c()
  
  for (plateN in 1:length(plateDataset)) {
    curr_plate <- c()
    curr_plate_mmScore <- c()
    tempData <- plateDataset[[plateN]]
    wellMean <- mean(tempData$TRG, na.rm =T)
    wellSD <- sd(tempData$TRG, na.rm = T)
    tempDataMin <- min(tempData$TRG, na.rm = T)
    tempDataMax <- max(tempData$TRG, na.rm = T)
    
    for (well in c(tempData$TRG)) {
      std_score <- (well - wellMean)/wellSD
      minMaxScore <- (well-tempDataMin)/(tempDataMax-tempDataMin)
      curr_plate <- c(curr_plate, std_score)
      curr_plate_mmScore <- c(curr_plate_mmScore, minMaxScore)
    }
    tempData$ZScore <- curr_plate
    tempData$mmScore <- curr_plate_mmScore
    std_rep <- c(std_rep, list(tempData))
  }
  
  
  return(std_rep)
}

dataEditor <- function(ExperimentPlate, Lib_Plate_Num, Plate_Treatments) {
  for (plateN in 1:length(ExperimentPlate)) {
    ExperimentPlate[[plateN]]$lpNum <- paste0('LP', Lib_Plate_Num)
    ExperimentPlate[[plateN]]$PlateTreatment <- c('UNT', Plate_Treatments, Plate_Treatments, Plate_Treatments)[plateN] 
    ExperimentPlate[[plateN]]$compoundID <- paste0(ExperimentPlate[[plateN]]$lpNum, '_', ExperimentPlate[[plateN]]$wellName)
    ExperimentPlate[[plateN]]$ConditionTreatment <- paste0(ExperimentPlate[[plateN]]$PlateTreatment, '_', ExperimentPlate[[plateN]]$lpNum, '_', ExperimentPlate[[plateN]]$wellName)
    
    ExperimentPlate[[plateN]] <- merge(x = ExperimentPlate[[plateN]], y = cleanLibMap, by = 'compoundID')
  }
  
  return(ExperimentPlate)
}

medNormFunc <- function(strainSub){
  dataMAD <- abs(mad(strainSub, na.rm = T))
  dataMED <- median(strainSub, na.rm = T)
  
  STDRES <- strainSub/dataMED
  return(STDRES)
}

lowessConverter <- function(dataset){
  tempDF <- data.frame('x' = dataset[[1]])
  finalDF <- data.frame('Time' = dataset[[1]])
  colIDs <- colnames(dataset)
  
  par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
  for (value in 2:length(dataset)) {
    tempDF[['y']] <- dataset[[value]]
    plot(tempDF, pch = 19, cex = 0.75, ylim= c(0,0.5), xlim = c(0,21))
    lowessFit <- lowess(tempDF, f = 0.6, iter=1)
    lowMyData <- data.frame('x' = lowessFit$x, 'y' = lowessFit$y)
    lines(lowMyData)
    
    finalDF[[colIDs[value]]] = lowessFit$y
  }
  #TRG_CALCULATOR_THREE(finalDF, date = NULL, plate_no = NULL)
  
  dev.off()
  return(finalDF)
}

serialPlotter <- function(whereIsFileLocated, 
                          saveFileAs, 
                          plateType,
                          sepStat = 'tog',
                          Treatment = 'Treatment',
                          Lib_Plate_Num,
                          untPos = 1,
                          smoothing = F) { 
  
  All_files = list.files(paste0(whereIsFileLocated, collapse = NULL))
  Plates = All_files[grep("txt",All_files)]
  n <- 0 
  my_curr_plates <- c()
  
  
  if (smoothing == F) {
    for (plate in Plates) {
      n <- n+1
      my_curr_ <- read.table(paste0(whereIsFileLocated, plate), header = T)
      my_curr_plates <- c(my_curr_plates, list(my_curr_))
    }
  } else {
    for (plate in Plates) {
      n <- n+1
      my_curr_ <- read.table(paste0(whereIsFileLocated, plate), header = T)
      my_curr_plates <- c(my_curr_plates, list(lowessConverter(my_curr_)))
    }
  }
  
  pdf(file = paste0(whereIsFileLocated, saveFileAs) ,useDingbats = F, paper = 'a4r', width = 120, height = 80)
  #Replicate Correlation
  par(mfrow = c(1,1))
  plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
  text(0,0, 'Replicate Correlation')
  
  n <- 1
  colIDs <- c('Plate1')
  melted_mtx <- data.frame('Plate1' = melt(my_curr_plates[[1]][2:length(my_curr_plates[[1]])])$value)
  for (currNo in 2:length(my_curr_plates)) {
    n <- n+1
    colIDs  <- c(colIDs, paste0('Plate', n))
    print(colIDs)
    curr <- my_curr_plates[[currNo]]
    newCol <- melt(curr[2:length(curr)])$value
    melted_mtx[colIDs[currNo]] <- newCol
  }
  
  pairs(melted_mtx, col = 'black', cex= 0.75, pch=21, bg='lightgreen', upper.panel = panel.cor, lower.panel = panel.smooth)
  
  ####create raw plots
  par(mfrow = c(1,1))
  plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
  text(0,0, 'Experiment 1')
  
  par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
  wellNames <- names(my_curr_plates[[1]])[2:length(my_curr_plates[[1]])]
  
  if (sepStat == 'tog') {
    for (well in 2:length(my_curr_plates[[1]])) {
      col_n <- 0
      plot(0,0, xlim = c(0,25), ylim = c(0,1.2), col='white', main = names(my_curr_plates[[1]])[well])
      legend('topleft', legend = c('unt', 'rep1','rep2','rep3'), col = 1:4, pch = 19, cex = 0.75)
      for (plate in my_curr_plates) {
        #print(names(plate[well]))
        col_n <- col_n+1
        points(x = plate[['Time']], y = plate[[well]], col = col_n)
      }
    }
  } else {
    for (plate in my_curr_plates) {
      col_n <- 0
      for (well in wellNames) {
        plot(plate[['Time']], plate[[well]], col = 4, pch = 19, cex = 0.75, 
             main = well, xlim = c(0,25), ylim = c(0,1.2))
      }
    }
  }
  
  
  
  # compute time of regrowth for each plate completely
  all_Data_HOLDER <- c()
  
  par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
  for (plate in my_curr_plates) {
    n = 0
    TEMP_DATA_HOLDER <- TRG_CALCULATOR_THREE(plate = plate, date = NULL, 
                                             plate_no = paste0('Plate ', n), 
                                             adjP = 0.15)
    
    all_Data_HOLDER <- c(all_Data_HOLDER, list(TEMP_DATA_HOLDER))
    n <- n+1
  }
  
  all_Data_HOLDER <- zScorer(all_Data_HOLDER)
  
  #take only TRG for plotting (now min-max normalized)
  TRG_HOLDER <- c()
  for (holding in all_Data_HOLDER) {
    TRG_HOLDER <- c(TRG_HOLDER, list(holding$TRG))
  }
  
  #take only TRG for plotting (now min-max normalized)
  NORM_TRG_HOLDER <- c()
  for (holding in all_Data_HOLDER) {
    NORM_TRG_HOLDER <- c(NORM_TRG_HOLDER, list(holding$mmScore))
  }
  
  #take only GRate for plotting
  GRATE_HOLDER <- c()
  for (holding in all_Data_HOLDER) {
    GRATE_HOLDER <- c(GRATE_HOLDER, list(holding$GRate))
  }
  #ignore for testing plates
  
  if (plateType == 384) {
    all_Data_HOLDER <- dataEditor(ExperimentPlate = all_Data_HOLDER, Lib_Plate_Num = Lib_Plate_Num, Plate_Treatments = Treatment)
  } else {
    all_Data_HOLDER <- all_Data_HOLDER
  }
  ###all_Data_HOLDER <- dataEditor(ExperimentPlate = all_Data_HOLDER, Lib_Plate_Num = Lib_Plate_Num, Plate_Treatments = Treatment)
  
  
  #PLOT PLATES WITH PLATEPLOTTER
  N <- 0
  for (compTRG in TRG_HOLDER) {
    N = N+1
    ppll <- wellPlatePlotter(SingleColData = c(compTRG), 
                             plateType = plateType, labelTitle = paste0('plate', N), colorScale = c(0,20))
    
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
      lims(x = c(0,25), y = c(0,2.5))
  }
  
  do.call(grid.arrange,NPLOTS)
  
  print(1)
  #outliers with replicates
  par(mfrow = c(round(length(TRG_HOLDER)/2), 2))
  HITS <- Reduce(intersect, list(Boxplot(all_Data_HOLDER[[2]]$TRG),
                                 Boxplot(all_Data_HOLDER[[3]]$TRG),
                                 Boxplot(all_Data_HOLDER[[4]]$TRG)))
  print(HITS)
  
  OUTLIER_DF <- c()
  for (plateOrder in 2:length(all_Data_HOLDER)) {
    OUTLIER_DF <- c(OUTLIER_DF , list(all_Data_HOLDER[[plateOrder]][HITS, ]))
  }
  
  ALL_OUTLIER_DATA_HOLDER <- c()
  for (plateOrder in 1:length(all_Data_HOLDER)) {
    ALL_OUTLIER_DATA_HOLDER <- c(ALL_OUTLIER_DATA_HOLDER, list(Boxplot(all_Data_HOLDER[[plateOrder]]$TRG)))
  }
  
  ALL_OUTLIER_DATA_LABEL_DF <- c()
  for (plateOrder in 1:length(all_Data_HOLDER)) {
    ALL_OUTLIER_DATA_LABEL_DF <- c(ALL_OUTLIER_DATA_LABEL_DF, list(all_Data_HOLDER[[plateOrder]][ALL_OUTLIER_DATA_HOLDER[[plateOrder]],]))
  }
  
  #Reduce(intersect, ALL_OUTLIER_DATA_HOLDER)
  
  OUTLIER_LABEL_HOLDER <- c()
  for (plateOrder in 1:length(OUTLIER_DF)) {
    OUTLIER_LABEL_HOLDER <- c(OUTLIER_LABEL_HOLDER, 
                              list(OUTLIER_DF[[plateOrder]]$wellName))
  }
  
  ALL_OUTLIER_LABEL_HOLDER <- c()
  for (plateOrder in 1:length(ALL_OUTLIER_DATA_LABEL_DF)) {
    ALL_OUTLIER_LABEL_HOLDER <- c(ALL_OUTLIER_LABEL_HOLDER, 
                                  list(ALL_OUTLIER_DATA_LABEL_DF[[plateOrder]]$wellName))
  }
  
  print(ALL_OUTLIER_LABEL_HOLDER)
  
  #PLOT VENN DIagram for intersection of outliers
  if (sepStat == 'tog') {
    veggven <- ggvenn(data = ALL_OUTLIER_LABEL_HOLDER[2:4], columns = 1:3, stroke_size = 0,
                      fill_color = c('#F6992F', '#24749B', '#22B258'), fill_alpha = 0.5, text_size = 5, show_elements = F)
    grid.arrange(veggven)
  }
  
  #boxplot to emphasize outliers
  print('completed ggven')
  #all_intersects <- intersect(intersect(OUTLIER_LABEL_HOLDER[[1]], OUTLIER_LABEL_HOLDER[[3]]), OUTLIER_LABEL_HOLDER[[4]])
  
  
  par(mfrow = c(1,1))
  plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
  text(0,0, toString(HITS))
  
  
  if (plateType == 384) {
    for (plate in 1:length(OUTLIER_LABEL_HOLDER)) {
      par(mfrow = c(1,1))
      plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
      mtext(paste0('Replicate ', plate))
      cTemPlate <- subset(all_Data_HOLDER[[plate+1]], wellName %in% OUTLIER_LABEL_HOLDER[[plate]])
      if (nrow(cTemPlate) > 0 ) {
        cTemPlate <- cTemPlate[c("Catalog.Number","TRG", "GRate", "Identifier", "ConditionTreatment", "ZScore")]
        grid.table(cTemPlate)
      }
      
    }
    
    for (plate in 1:length(ALL_OUTLIER_LABEL_HOLDER)) {
      par(mfrow = c(1,1))
      plot(0,0, col = 0, xaxt='n', yaxt='n', ann=FALSE)
      mtext(paste0('All Outliers: Replicate ', plate))
      cTemPlate <- subset(all_Data_HOLDER[[plate]], wellName %in% ALL_OUTLIER_LABEL_HOLDER[[plate]])
      if (nrow(cTemPlate) > 0 ) {
        cTemPlate <- cTemPlate[c("Catalog.Number","TRG", "GRate", "Identifier", "ConditionTreatment", "ZScore")]
        grid.table(cTemPlate)
      }
    }
    
    
    if (sepStat == 'tog') {
      par(mfrow = c(1,1), mar = c(5, 5, 5, 5))
      
      #boxplot of all data
      bxpdat <- boxplot(TRG_HOLDER, ylim= c(0,25), outline = F, varwidth = T, col = 'grey', names = paste0('Plate ', 1:length(TRG_HOLDER)),
                        boxwex=.8, notch = T, border = 1, lwd=2, cex.lab = 1.5, cex.axis = 1.2,
                        ylab = substitute(paste(bold('Time of Regrowth (Hrs)'))), 
                        xlab=substitute(paste(bold(Treatment))))
      print('completed box plots')
      
      
      n <- 0
      for (labelling in 1:length(OUTLIER_LABEL_HOLDER)) {
        bxpdat <- boxplot(TRG_HOLDER[[1]], plot = F)
        n <- n+1
        text(x = 1+labelling, y= bxpdat$out, labels = c(OUTLIER_LABEL_HOLDER[[labelling]], ''), col = 1:10, cex=0.75) 
      }
    }
    
    #SCORING
    for (plateN in 1:length(all_Data_HOLDER)) {
      TRG_NULL <- mean(subset(all_Data_HOLDER[[plateN]], Catalog.Number == 'DMSO')$TRG, na.rm = T)
      print(TRG_NULL)
      TRG_DIFF <- c() #wellTRG -dmsoWellTrg (trgNull)
      TRG_DIV <- c() #wellTRG / dmsoWellTrg (trgNull)
      for (wellN in 1:length(all_Data_HOLDER[[plateN]]$TRG)) {
        TRG_DIFF <- c(TRG_DIFF, all_Data_HOLDER[[plateN]]$TRG[wellN] - TRG_NULL)
        TRG_DIV <- c(TRG_DIV, all_Data_HOLDER[[plateN]]$TRG[wellN] / TRG_NULL)
      }
      all_Data_HOLDER[[plateN]]$TRG_DIFF <- TRG_DIFF
      all_Data_HOLDER[[plateN]]$TRG_DIV <- TRG_DIV
    }
    
    UntreatedPlate <- all_Data_HOLDER[[untPos]] #hardcoded... not the best, but will fix, the last plates
    UntreatedPlateDiv <- all_Data_HOLDER[[untPos]] 
    
    for (plateN in 1:length(all_Data_HOLDER)) {
      all_Data_HOLDER[[plateN]]$DIFF_SPACER <- all_Data_HOLDER[[plateN]]$TRG_DIFF - UntreatedPlate$TRG_DIFF
      all_Data_HOLDER[[plateN]]$DIV_SPACER <- all_Data_HOLDER[[plateN]]$TRG_DIV / UntreatedPlateDiv$TRG_DIV
    }
    
    #median adjusted
    for (plateN in 1:length(all_Data_HOLDER)) {
      all_Data_HOLDER[[plateN]]$medAdjTRG <- with(all_Data_HOLDER[[plateN]], medNormFunc(TRG))
    }
    
    for (plateN in 1:length(all_Data_HOLDER)) {
      med_TRG_NULL <- mean(subset(all_Data_HOLDER[[plateN]], Catalog.Number == 'DMSO')$medAdjTRG, na.rm = T)
      print(med_TRG_NULL)
      med_TRG_DIFF <- c() #wellTRG -dmsoWellTrg (trgNull)
      med_TRG_DIV <- c() #wellTRG / dmsoWellTrg (trgNull)
      for (wellN in 1:length(all_Data_HOLDER[[plateN]]$medAdjTRG)) {
        med_TRG_DIFF <- c(med_TRG_DIFF, all_Data_HOLDER[[plateN]]$medAdjTRG[wellN] - med_TRG_NULL)
        med_TRG_DIV <- c(med_TRG_DIV, all_Data_HOLDER[[plateN]]$medAdjTRG[wellN] / med_TRG_NULL)
      }
      all_Data_HOLDER[[plateN]]$med_TRG_DIFF <- med_TRG_DIFF
      all_Data_HOLDER[[plateN]]$med_TRG_DIV <- med_TRG_DIV
    }
    
    UntreatedPlate <- all_Data_HOLDER[[untPos]] #hardcoded... not the best, but will fix, the last plates
    UntreatedPlateDiv <- all_Data_HOLDER[[untPos]] 
    
    for (plateN in 1:length(all_Data_HOLDER)) {
      all_Data_HOLDER[[plateN]]$med_DIFF_SPACER <- all_Data_HOLDER[[plateN]]$med_TRG_DIFF - UntreatedPlate$med_TRG_DIFF
      all_Data_HOLDER[[plateN]]$med_DIV_SPACER <- all_Data_HOLDER[[plateN]]$med_TRG_DIV / UntreatedPlateDiv$med_TRG_DIV
    }
  } else {
    all_Data_HOLDER <- all_Data_HOLDER
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


#######---- Load data ------
setwd("~/Documents/wale_docs/phd/data/")
load('Phase_1_Data/Full_Data/Amoxicilin_Data.RData')
amoxData <- amikData
load('Phase_1_Data/Full_Data/Ciprofloxacin_Data.RData')
ciproData <- amikData
load('Phase_1_Data/Full_Data/Amikacin_Data.RData')

#######---- Hit detection scripts ------
distiller <- function(fullDataset, ColToSelect) {
  plateOrder <- c('unt','rep1','rep2','rep3')
  tempdf <- NULL
  for (lp in fullDataset) {
    
    n <- 1
    for (plate in lp) {
      newGen <- plate[c(ColToSelect, 'compoundID', 'lpNum', 'PlateTreatment','wellName')]
      newGen$UNT <- lp[[1]][[ColToSelect]]
      newGen$ORDER <- plateOrder[n]
      tempdf <- rbind(tempdf, newGen)
      n <- n+1
    }
  }
  return(tempdf)
}

distiller4plots <- function(fullDataset, ColToSelect) {
  plateOrder <- c('unt','rep1','rep2','rep3')
  tempdf <- NULL
  for (lp in fullDataset) {
    
    n <- 1
    for (plate in lp) {
      newGen <- plate[c(ColToSelect, 'compoundID', 'lpNum', 'PlateTreatment','wellName')]
      newGen$UNT <- lp[[1]][[ColToSelect]]
      newGen$ORDER <- plateOrder[n]
      tempdf <- rbind(tempdf, newGen)
      n <- n+1
    }
  }
  tempdf <- subset(tempdf, ORDER != 'unt')
  return(tempdf)
}

distIpper <- function(fullDataset, ColToSelect) {
  plateOrder <- c('unt','rep1','rep2','rep3')
  tempdf <- NULL
  for (lp in fullDataset) {
    
    
    
    newGen <- lp[[1]][c('compoundID', 'lpNum', 'wellName')]
    newGen$UNT <- lp[[1]][[ColToSelect]]
    newGen$REP1 <- lp[[2]][[ColToSelect]]
    newGen$REP2 <- lp[[3]][[ColToSelect]]
    newGen$REP3 <- lp[[4]][[ColToSelect]]
    tempdf <- rbind(tempdf, newGen)
    
    
  }
  return(tempdf)
}

trueOutlierDetectives <- function(fullDataset, n = 2) {
  
  untcut <- 1.5*IQR(fullDataset$UNT, na.rm = T) + 
    quantile(fullDataset$UNT, na.rm = T)[[4]] 
  
  rep1cut <- (1.5*IQR(fullDataset$REP1, na.rm = T) + 
                quantile(fullDataset$REP1, na.rm = T)[[4]])*n 
  
  rep2cut <- (1.5*IQR(fullDataset$REP2, na.rm = T) + 
                quantile(fullDataset$REP2, na.rm = T)[[4]])*n 
  
  rep3cut <- n*(1.5*IQR(fullDataset$REP3, na.rm = T) + 
                  quantile(fullDataset$REP3, na.rm = T)[[4]]) 
  
  
  untOut <- subset(fullDataset, UNT > untcut)$compoundID
  rep1Out <- subset(fullDataset, REP1 > rep1cut)$compoundID
  rep2Out <- subset(fullDataset, REP2 > rep2cut)$compoundID
  rep3Out <- subset(fullDataset, REP3 > rep3cut)$compoundID
  
  cut1 <- setdiff(rep1Out, untOut)
  cut2 <- setdiff(rep2Out, untOut)
  cut3 <- setdiff(rep3Out, untOut)
  
  opList = c(subset(data.frame(table(c(cut1, cut2, cut3)), stringsAsFactors = F), Freq >= 2)$Var1)
  
  return(as.vector(opList))
}

####### DMSO correlation #####
corrDataTableDMSO <- function(whereIsFileLocated) {
  
  All_files = list.files(paste0(whereIsFileLocated, collapse = NULL))
  Plates = All_files[grep("txt",All_files)]
  
  #import plates
  my_curr_plates <- c()
  
  n <- 0
  for (plate in Plates) {
    n <- n+1
    my_curr_ <- read.table(paste0(whereIsFileLocated, plate), header = T)[allDmsoWell]
    my_curr_plates <- c(my_curr_plates, list(my_curr_))
  }
  
  n <- 1
  colIDs <- c('Plate1')
  melted_mtx <- data.frame('Plate1' = melt(my_curr_plates[[1]][2:length(my_curr_plates[[1]])])$value)
  for (currNo in 2:length(my_curr_plates)) {
    n <- n+1
    colIDs  <- c(colIDs, paste0('Plate', n))
    print(colIDs)
    curr <- my_curr_plates[[currNo]]
    newCol <- melt(curr[2:length(curr)])$value
    melted_mtx[colIDs[currNo]] <- newCol
  }
  repMelts <- melted_mtx[1:4]
  #select just waht you need here
  whatIneed <- c(cor(repMelts))
  return(whatIneed)
}


#function for replicate correlations
corrDataTable <- function(whereIsFileLocated) {
  
  All_files = list.files(paste0(whereIsFileLocated, collapse = NULL))
  Plates = All_files[grep("txt",All_files)]
  
  #import plates
  my_curr_plates <- c()
  
  for (plate in Plates) {
    n <- n+1
    my_curr_ <- read.table(paste0(whereIsFileLocated, plate), header = T)
    my_curr_plates <- c(my_curr_plates, list(my_curr_))
  }
  
  n <- 1
  colIDs <- c('Plate1')
  melted_mtx <- data.frame('Plate1' = melt(my_curr_plates[[1]][2:length(my_curr_plates[[1]])])$value)
  for (currNo in 2:length(my_curr_plates)) {
    n <- n+1
    colIDs  <- c(colIDs, paste0('Plate', n))
    print(colIDs)
    curr <- my_curr_plates[[currNo]]
    newCol <- melt(curr[2:length(curr)])$value
    melted_mtx[colIDs[currNo]] <- newCol
  }
  repMelts <- melted_mtx[2:4]
  #select just waht you need here
  whatIneed <- c(cor(repMelts))[c(2,3,6)]
  return(whatIneed)
}

########
#plot across replicates at once

######------ plot replicates across a lp at once ######------ 
replicatePlotter <- function(whereIsFileLocated, well) {
  
  All_files = list.files(paste0(whereIsFileLocated, collapse = NULL))
  Plates = All_files[grep("txt",All_files)]
  
  myFilesLocated <- list()
  for (file in Plates) {
    myFilesLocated <- c(myFilesLocated, list(read.table(paste0(whereIsFileLocated, file), header = T)))
  } 
  
  plot(x = myFilesLocated[[1]]$Time, y = myFilesLocated[[1]][[well]], 
       xlim = c(1,24), ylim = c(0,1.3), type = 'l', xlab = 'Time', ylab = 'OD600',
       lwd = 4)#, main = paste0('Compound ', well))
  
  #plot DMSO4 (A1) untreated
  lines(x = myFilesLocated[[1]]$Time, y = myFilesLocated[[1]][['D21']],
        lwd = 4, col = 'grey')
  
  #plot DMSO4 (A1) treated
  lines(x = myFilesLocated[[1]]$Time, y = myFilesLocated[[3]][['D21']],
        lwd = 4, col = 'chocolate')
  
  #plot the wells in replicates
  for (repPla in 2:length(myFilesLocated)) {
    lines(x = myFilesLocated[[repPla]]$Time, y= myFilesLocated[[repPla]][[well]],
          lwd = 4, col = repPla)
  }
  
  #plot the legend
  legend('bottomright', legend = c(paste0('Untreated+',well), paste0(well, ' rep1'), 
                                   paste0(well, ' rep2'), paste0(well, ' rep3'),
                                   paste0('Untreated+DMSO'), paste0('Treated+DMSO')), col = c(1:4, 'grey', 'chocolate'), 
         pch = 19, cex = 0.5)
  
}

######------ modelling with david's model ######------ 
odDetection <- 0.01113222
odDetectoR <- function(fullDataset) {
  odDetection <- 0.01113222
  tempData <- list()
  for (plate in fullDataset) {
    plate$X_INTERCEPT <- (-1*plate$Y_INTERCEPT)/plate$GRate
    plate$detectionIntercept <- (odDetection - plate$Y_INTERCEPT)/plate$GRate
    tempData <- c(tempData, list(plate))
  }
  return(tempData)
}
#apply to all
amikData = lapply(amikData, odDetectoR)
amoxData = lapply(amoxData, odDetectoR)
ciproData = lapply(ciproData, odDetectoR)

diffModeller <- function(fullDataset) {
  tempData <- list()
  for (plate in fullDataset) {
    plate$X_intercept_diff <- plate$detectionIntercept - mean(subset(plate, wellName %in% allDmsoWell)$detectionIntercept, na.rm =T)
    plate$Y_intercept_diff <- plate$Y_INTERCEPT - mean(subset(plate, wellName %in% allDmsoWell)$Y_INTERCEPT, na.rm = T)
    tempData <- c(tempData, list(plate))
  }
  return(tempData)
}

amikData = lapply(amikData, diffModeller)
amoxData = lapply(amoxData, diffModeller)
ciproData = lapply(ciproData, diffModeller)


#predict
getXiNtercept<- function(fullDataset) {
  tempData <- list()
  
  unt <- fullDataset[[1]]$X_intercept_diff #think about this
  unty <- fullDataset[[1]]$Y_intercept_diff
  
  for (plate in fullDataset) {
    plate$PredXInterceptTreated <- unt+ plate$detectionIntercept
    plate$PredYInterceptTreated <- unty+ plate$detectionIntercept
    tempData <- c(tempData, list(plate))
  }
  return(tempData)
}


amikData = lapply(amikData, getXiNtercept)
amoxData = lapply(amoxData, getXiNtercept)
ciproData = lapply(ciproData, getXiNtercept)


#using david's formula
dvForm <- function(fullDataset) {
  tempData <- list()
  odDetection <- 0.01113222
  unt <- fullDataset[[1]]
  for (plate  in fullDataset) {
    plate$davForm <- (log(odDetection)-(subset(plate, wellName == 'A1')$Y_INTERCEPT + unt$Y_INTERCEPT - subset(unt, wellName == 'A1')$Y_INTERCEPT))/plate$GRate
    
    tempData <- c(tempData, list(plate))
  }
  return(tempData)
}

amikData = lapply(amikData, dvForm)
amoxData = lapply(amoxData, dvForm)
ciproData = lapply(ciproData, dvForm)

#plot the schema
plotIntercepts <- function(lp, repl = 1) {
  plateToPick = repl+1
  
  for (rowID in 1:384) {
    plot(0,0, xlim= c(0,20), ylim = c(-20,10), col = 'white')
    with(subset(lp[[1]], wellName == 'A1'), abline(a = Y_INTERCEPT, b = GRate))
    if (!is.na(lp[[1]][rowID,]$TRG)) {
      with(lp[[1]][rowID,], abline(a = Y_INTERCEPT, b = GRate))
    } else {
      abline(0,0)
    }
    
    with(subset(lp[[plateToPick]], wellName == 'A1'), abline(a = Y_INTERCEPT, b = GRate, col = plateToPick))
    if (!is.na(lp[[plateToPick]][rowID,]$TRG)) {
      with(lp[[plateToPick]][rowID,], abline(a = Y_INTERCEPT, b = GRate, col = plateToPick))
    } else {
      abline(0,0)
    }
    
    abline(v =0)
    abline(h =odDetection)
  }
}

##
######## All well names #####
allWellNames <- names(read.table(file = 'Phase_1_Data/Full_Data/Full_Amikacin/lp1/Plate1.txt', header = T))

#### DMSO SELECTOR #####

dmsoSelector <- function(fullDataset, treatment_No, ColToSelect) {
  #1 = unt
  #2,3,4 = rep1, rep2, rep3
  
  tempDF <- unique(fullDataset[[1]][[treatment_No]])['wellName']
  tempDF$LP1 <- unique(fullDataset[[1]][[treatment_No]])[ColToSelect]
  tempDF$LP2 <- fullDataset[[2]][[treatment_No]][ColToSelect]
  tempDF$LP3 <- fullDataset[[3]][[treatment_No]][ColToSelect]
  tempDF$LP4 <- fullDataset[[4]][[treatment_No]][ColToSelect]
  tempDF$LP5 <- fullDataset[[5]][[treatment_No]][ColToSelect]
  tempDF$LP6 <- fullDataset[[6]][[treatment_No]][ColToSelect]
  
  tempDF <- subset(tempDF, wellName %in% allDmsoWell)
  
  return(tempDF)
}

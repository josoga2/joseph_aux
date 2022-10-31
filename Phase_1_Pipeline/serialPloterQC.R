#accepts growth curve data (od over time) in txt format in the same folder
#returns a qc report of the dataset with preliminary analysis
#whereIsFileLocated: working directory of the files (string)
#saveFileAs: output pdf file name of qc report (string)
#plateType: 96/384 (numeric)
#sepstat: 'tog'-together; comapre the files together 
#Treatment: antibiotic used or condition (string)
#Lib_Plate_Num: number of the library plate used: 1 (int), can be ignored as NULL
#untPos: location of the untreated/control dataset file
#smoothing: whether to smoothen or not with lowess
#author: Adewale Joseph Ogunleye


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
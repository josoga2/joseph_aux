#plateDataset: a list of 384 or 96 well plates of any length
#ColToSelect: name of the column to calculate z score on
#adds 2 column to plateDataset: [1. Zscore: x-mean(x)/std_dev(x) ], [2: x-min(x)/max(x)-min(x)]
#author: Adewale Joseph Ogunleye


zScorer <- function(plateDataset, ColToSelect) {
  std_rep <- c()
  
  for (plateN in 1:length(plateDataset)) {
    curr_plate <- c()
    curr_plate_mmScore <- c()
    tempData <- plateDataset[[plateN]]
    wellMean <- mean(tempData[[ColToSelect]], na.rm =T)
    wellSD <- sd(tempData[[ColToSelect]], na.rm = T)
    tempDataMin <- min(tempData[[ColToSelect]], na.rm = T)
    tempDataMax <- max(tempData[[ColToSelect]], na.rm = T)
    
    for (well in c(tempData[[ColToSelect]])) {
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

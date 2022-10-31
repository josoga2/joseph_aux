#set of functions for modeling TRG between treatment and untreated datasets
#author: Adewale Joseph Ogunleye

odDetection <- 0.01113222 #tunable
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
#amikData = lapply(amikData, odDetectoR)
#amoxData = lapply(amoxData, odDetectoR)
#ciproData = lapply(ciproData, odDetectoR)

diffModeller <- function(fullDataset) {
  tempData <- list()
  for (plate in fullDataset) {
    plate$X_intercept_diff <- plate$detectionIntercept - mean(subset(plate, wellName %in% allDmsoWell)$detectionIntercept, na.rm =T)
    plate$Y_intercept_diff <- plate$Y_INTERCEPT - mean(subset(plate, wellName %in% allDmsoWell)$Y_INTERCEPT, na.rm = T)
    tempData <- c(tempData, list(plate))
  }
  return(tempData)
}

#amikData = lapply(amikData, diffModeller)
#amoxData = lapply(amoxData, diffModeller)
#ciproData = lapply(ciproData, diffModeller)


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


#amikData = lapply(amikData, getXiNtercept)
#amoxData = lapply(amoxData, getXiNtercept)
#ciproData = lapply(ciproData, getXiNtercept)


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
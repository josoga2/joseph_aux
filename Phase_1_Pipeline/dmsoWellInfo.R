
#select only DMSO wells
#author: Adewale Joseph Ogunleye

######## DMSO4 DATA #####
#pick all DMSO wells
LibMap <- read.table('../PLANS/LibMap.txt', sep = '\t', quote = '\t', header = T)
allDmsoWell <- unique(subset(LibMap, Catalog.Number == 'DMSO')$Well)
#########################

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
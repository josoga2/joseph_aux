#fullDataset: dataset containing a list of datasets arranged by library plates. each lp is arranged as a list of treamtents
#Coltoselect: the name of desired column/feature to be exported/sliced
#returns: a distilled/melted version of dataset sorted by library plate number and treatments (untreated is duplicated)
#author: Adewale Joseph Ogunleye


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
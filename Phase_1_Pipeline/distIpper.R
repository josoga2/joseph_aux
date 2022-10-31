#fullDataset: dataset containing a list of datasets arranged by library plates. each lp is arranged as a list of treamtents
#Coltoselect: the name of desired column/feature to be exported/sliced
#returns: a distilled/melted version of dataset with seprate columns for each replicate, but sorted by library plate number and treatments 
#author: Adewale Joseph Ogunleye

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
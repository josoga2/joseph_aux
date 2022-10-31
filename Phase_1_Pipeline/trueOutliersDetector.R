#fullDataset: accepts distipper output
#n: size factor of intercepts calculation
#returns: a vector of outlier row identifiers
#author: Adewale Joseph Ogunleye

trueOutlierDetectives <- function(fullDataset, n = 1.5) {
  
  untcut <- 1.5*IQR(fullDataset$UNT, na.rm = T) + 
    quantile(fullDataset$UNT, na.rm = T)[[4]] 
  
  rep1cut <- n*(1.5*IQR(fullDataset$REP1, na.rm = T) + 
                quantile(fullDataset$REP1, na.rm = T)[[4]])*n 
  
  rep2cut <- n*(1.5*IQR(fullDataset$REP2, na.rm = T) + 
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
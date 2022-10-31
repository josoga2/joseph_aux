#for lowess smoothing 
# can be edited to return plots and trg calculation result
#accepts any dataframe with first column as Time (not hard rule, but needed for visualization)
#returns lowess smoothened result 
#author: Adewale Joseph Ogunleye


lowessConverter <- function(dataset){
  tempDF <- data.frame('x' = dataset[[1]])
  finalDF <- data.frame('Time' = dataset[[1]])
  colIDs <- colnames(dataset)
  
  #par(mfrow = c(8,12), cex=0.25, mar=c(2,2,2,2), oma=c(2,2,2,2), no.readonly = T)
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
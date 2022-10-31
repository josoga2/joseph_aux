#calculates median normalization
#author: Adewale Joseph Ogunleye

medNormFunc <- function(strainSub){
  dataMAD <- abs(mad(strainSub, na.rm = T))
  dataMED <- median(strainSub, na.rm = T)
  
  STDRES <- strainSub/dataMED
  return(STDRES)
}
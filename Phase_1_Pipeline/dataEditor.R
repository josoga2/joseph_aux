#for annotating TRG_CALCULATOR OUTPUT with library plate order
#accepts 
#ExperimentPlate: a list of 384 or 96 well plates of length 4
#Lib_Plate_Num: number of library plate (numeric)
#Plate_Treatments: string, usually untreated, antibiotics, dmso, etc

#adds 4 columns to ExperimentPlate: 
#1. concat str of 'LP'+Lib_Plate_Num, 
#2. Plate_Treatments, 
#3. concat str of Lib_Plate_Num+_+wellName
#4. condition
#author: Adewale Joseph Ogunleye

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
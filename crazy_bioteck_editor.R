#edit biotek output
#for only crazy *1000 outputs

editor_1000 <- function(data_location, editedFileName) {
  df <- read.table(file = data_location, sep = '\t', header = T)
  editedFileName <- paste0(data_location, editedFileName)
  df <- cbind(df[1],df[2:ncol(df)]/1000)
  write.table(df, file = editedFileName, sep = '\t', row.names = F, quote = F)
  return(df)
}

#implementation
wer = editor_1000(data_location = '/Users/josoga2/Documents/wale_docs/phd/data/Phase_1_Data/23_07_22/Output/23_07_22_single_pin.txt', editedFileName = 'edit')
head(wer)

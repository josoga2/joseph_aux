#visualizing fits using slope (grate) and yintercet information
#author: Adewale Joseph Ogunleye


plotIntercepts <- function(lp, repl = 1) {
  plateToPick = repl+1
  
  for (rowID in 1:384) {
    plot(0,0, xlim= c(0,20), ylim = c(-20,10), col = 'white')
    with(subset(lp[[1]], wellName == 'A1'), abline(a = Y_INTERCEPT, b = GRate))
    if (!is.na(lp[[1]][rowID,]$TRG)) {
      with(lp[[1]][rowID,], abline(a = Y_INTERCEPT, b = GRate))
    } else {
      abline(0,0)
    }
    
    with(subset(lp[[plateToPick]], wellName == 'A1'), abline(a = Y_INTERCEPT, b = GRate, col = plateToPick))
    if (!is.na(lp[[plateToPick]][rowID,]$TRG)) {
      with(lp[[plateToPick]][rowID,], abline(a = Y_INTERCEPT, b = GRate, col = plateToPick))
    } else {
      abline(0,0)
    }
    
    abline(v =0)
    abline(h =odDetection)
  }
}
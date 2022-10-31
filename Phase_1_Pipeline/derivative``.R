#second derivative calculator for unsmoothened growth curve (od over time)
#returns the second derivative of y (OD) over the same x(Time)
#input: a 2 column dataset, with one of the columns as wellName (df)
#pl: wether to plot second derviative result (boolean)
#adjparam: parameter for adjusting baseline search (numeric)

derivative <- function(input, pl=T, adjParam = 0.1){#you can set pl = False to reduce the number of plots generated in the main function
  
  minData <- min(input$wellName, na.rm = T)
  #print(minData)
  
  input <- subset(input, wellName < minData+adjParam)
  #print(length(input$Time))
  #print(input$wellName)
  grads= c()
  slope <- runner(
    x = input,
    k = 4,
    idx =  c(1:length(input$Time)),
    function(x) {
      coefficients(lm(wellName ~ Time, data = x))[[2]]
    }
  )
  
  grads <- c(grads, slope)
  #print(grads)
  
  
  roundGrads <- round(grads, digit=3)
  #print(roundGrads)
  maxGrads = NULL
  
  if (length(roundGrads) > 3) {
    selection <- roundGrads[4:length(roundGrads)]
    maxGrads = which(roundGrads == max(selection, na.rm = T))-4
  }else{
    maxGrads = which(roundGrads == max(roundGrads, na.rm = T))-4
  }
  
  
  
  secondInput <- data.frame('Time' = input$Time, 'grads'=grads)
  #print()
  
  secondGrads= c()
  slope2 <- runner(
    x = na.omit(secondInput),
    k = 4,
    #idx =  c(1:length(secondInput$grads)),
    function(x) {
      coefficients(lm(grads ~ Time, data = x))[[2]]
    }
  )
  
  secondGrads <- c(secondGrads, slope2)
  #print(secondGrads)
  #6
  #print(tail(secondGrads, 35))
  
  if (length(subset(input, wellName > minData, wellName < minData+0.2)$Time) > 6) {
    solution <- which(secondGrads == max(tail(secondGrads[5:(length(secondGrads)-10)], -10), na.rm = T))
  } else {
    solution <- which(secondGrads == max(tail(secondGrads, 6), na.rm = T))
  }
  
  
  OD <- NULL
  
  #print(solution)
  if (is.finite(solution) && length(solution) == 1) {
    OD <- input$wellName[solution]
    if (pl==T) {
      #plot(input$wellName, ylim = c(-0.2,0.4))
      #lines(secondGrads, col='blue')
      #abline(v = solution)
      #abline(h= OD)
    }
  } else {
    solution = 0
    OD = 0
  }
  
  
  return(c(solution, OD,maxGrads))
}
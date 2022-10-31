#accepts a variable that correspond to an integer or float value for wells in 
#a 384 or 96 well plate
#SingleColData: variable to be plotted over 384 or 96 wp
#plateType: 96 or 384 wp
#labelTitle: what to label the plot as
#colorScale (numeric): min and max of color scale to use for plotting 
#returns a plate-like plot of the data arranged as 8*12 or 16*24 
#author: Adewale Joseph Ogunleye

wellPlatePlotter <- function(SingleColData, plateType, labelTitle, colorScale = c(0,25)) {
  
  if (plateType == 96) {
    mymeta <- metafile96
    mymeta$type <- 'S1'
    mymeta$id <- 'Sample'
    mymeta$concentration <- 0.15
    rawdat <- cbind('X' = rawdata96$X , t(data.frame(matrix(data = SingleColData, nrow = 12, ncol = 8) )))
    #print(rawdat)
    rawdat <- data2plateformat(data = rawdat, platetype = 96)
    trg_df<- plate2df(rawdat)
    
    pll <- ggplot(data = trg_df, aes(x = reorder(x = row, desc(row)), y =col, col=value))+
      geom_point(size=18)+ylim(c(0,13))+coord_flip()+theme_bw() +
      scale_y_continuous(position = 'right' ,breaks=c(1:12), labels=c(1:12))+ 
      scale_color_gradient2(low = '#2a6a99', mid = '#c6dbef', midpoint = mean(colorScale), high = '#d88546', limits= colorScale)+
      labs(x=NULL, y = NULL)+ theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_text(label= round(c(trg_df$value), 2), color = 'black', size=5)+
      ggtitle(label = paste0('Plate '))
  }
  
  if (plateType == 384) {
    mymeta <- metafile384
    mymeta$type <- 'S1'
    mymeta$id <- 'Sample'
    mymeta$concentration <- 0.15
    rawdat <- cbind('X' = rawdata384$X , t(data.frame(matrix(data = SingleColData, nrow = 24, ncol = 16) )))
    #print(rawdat)
    rawdat <- data2plateformat(data = rawdat, platetype = 384)
    trg_df<- plate2df(rawdat)
    
    pll <- ggplot(data = trg_df, aes(x = reorder(x = row, desc(row)), y =col, col=value))+
      geom_point(size=12.5)+ylim(c(0,13))+coord_flip()+theme_bw() +
      scale_y_continuous(position = 'right' ,breaks=c(1:24), labels=c(1:24))+ 
      scale_color_gradient2(low = '#2a6a99', mid = '#c6dbef', midpoint = mean(colorScale), high = '#d88546', limits= colorScale)+
      labs(x=NULL, y = NULL)+ theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_text(label= round(c(trg_df$value), 2), color = 'black', size=3)+
      ggtitle(label = paste0( labelTitle))
  }
  return(pll)
}

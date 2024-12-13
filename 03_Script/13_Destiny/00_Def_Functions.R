
# Functions



exp2col<-function(value,breaks=1000,
                  cols=viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1)){
  colPal <- colorRampPalette(cols)
  marker.color <- colPal(breaks)[as.numeric(cut(value, breaks = breaks))]
  return(marker.color)
}


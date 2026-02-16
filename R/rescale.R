rescale<-function(x, a, b) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  a + (x - x_min) * (b - a) / (x_max - x_min)}


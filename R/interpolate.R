interpolate <-
function (obj, grid.list) {
  # Performs a linear interpolation of values in "obj" 
  # to a new series of points assigned by values in "grid"
  # Procedure is parallelised
  #Setup Framework
  x <- grid.list$x
  y <- grid.list$y
  M <- length(x)
  N <- length(y)
  #Initialise output matrix 
  out <- matrix(NA, nrow = M, ncol = N)
  # Perform interpolation
  out<-foreach(i=1:M, .inorder=TRUE) %dopar% {
    interp2D(obj, cbind(rep(x[i], N), y))
  }
  out<-array(unlist(out), dim=c(dim(out[[1]]),length(out)))
  #Return list of grid points and interpolated values
  return(as.array(list(x = x, y = y, z = out)))
}

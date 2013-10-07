makethetarho <-
function(dim) {
  # Produces a list containing a unit radial map and an angular map
  # of a 2D array. I.e. the radial map has values increasing with 
  # distance from the centre of the array; the angular map has elements 
  # equal to the angle from "north"
  message('-----------------------Make_Theta_Rho----------------------------------')

  #Set dimensions
  if (length(dim) >= 1) {
    w = dim[1]
    if (length(dim) == 2) {
      h = dim[2] 
    } else {
      h = w
    }
  } else { sink(type="message") ; stop("In makerho() - Dimensions not specified") }

  #Setup Grids
  x = (t(array(0:(w*h), dim=c(w,h))) %% w)-(w-1.)/2.
  y = (t(array(0:(w*h), dim=c(w,h)))/w)-(h-1)/2.
  #Perform Calcultions
  rho = sqrt(x^2+y^2)
  message('==========END==========Make_Theta_Rho==============END=================\n')
  theta = atan2(y,x)
  return(as.array(list(rho=rho,theta=theta)))
}

objpos2objpix <-
function(x_g, y_g) {
  #Find and return the Pixel (centre; x_p,y_p) closest to Gama position (x_y,y_g)
  message('--------------------------ObjPos_2_ObjPix------------------------------')
  #x_p<-round(x_g, digits=0)
  #y_p<-round(y_g, digits=0)
  x_p<-floor(x_g)
  y_p<-floor(y_g)

  #-----Diagnostic-----#
  if (diagnostic) { 
    message("X_g Y_g")
    sink(file=sinkfile, type="output", append=TRUE)
    print(summary(x_g))
    print(summary(y_g))
    message("X_p Y_p")
    print(summary(x_p))
    print(summary(y_p))
    sink(file=NULL, type="output")
  }
  message('===========END==========ObjPos_2_ObjPix============END=================\n')
  return(rbind(x_p, y_p))
}

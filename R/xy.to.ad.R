xy.to.ad <-
function(x,y,astr.struc,diagnostic=FALSE){

  if (!missing(x)) {
    if (missing(y) & length(dim(x))==2) {
      if (ncol(x)!=2) { stop(paste0('x must be of dim 2, not ',ncol(x))) }
      y<-x[,2]
      x<-x[,1]
    }
  }
  if (missing(x)) {
    x<-c(1,astr.struc$NAXIS[1])+0.5
  }
  if (missing(y)) {
    y<-c(1,astr.struc$NAXIS[2])+0.5
  }
  ra0=astr.struc$CRVAL[1]
  dec0=astr.struc$CRVAL[2]
  x0=astr.struc$CRPIX[1]
  y0=astr.struc$CRPIX[2]
  xscale=astr.struc$CD[1,1]
  yscale=astr.struc$CD[2,2]
  if (diagnostic) { message(paste("ASTR:",ra0,dec0,x0,y0,xscale,yscale)) }
  radec<-xy.to.radec(x,y,ra0,dec0,x0,y0,xscale,yscale)
  return=radec
}

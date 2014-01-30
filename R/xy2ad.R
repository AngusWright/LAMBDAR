xy2ad <-
function(x,y,astr_struc,diagnostic=FALSE){

  ra0=astr_struc$CRVAL[1]
  dec0=astr_struc$CRVAL[2]
  x0=astr_struc$CRPIX[1]
  y0=astr_struc$CRPIX[2]
  xscale=astr_struc$CD[1,1]
  yscale=astr_struc$CD[2,2]
  if (diagnostic) { message(paste("ASTR:",ra0,dec0,x0,y0,xscale,yscale)) }
  radec<-xy2radec(x,y,ra0,dec0,x0,y0,xscale,yscale)
  return=radec
}

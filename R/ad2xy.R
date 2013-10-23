ad2xy <-
function(ra,dec,astr_struc,diagnostic=FALSE){
# Converts RA/DEC (degrees) to x/y (pixels) position using the TAN Gnomonic projection system
# point of interest RA & DEC, anchor point RA & DEC, anchor point x & y, x & y scale (degrees per pixel)
# FROM: http://mathworld.wolfram.com/GnomonicProjection.html
  ra0=astr_struc$CRVAL[1]
  dec0=astr_struc$CRVAL[2]
  x0=astr_struc$CRPIX[1]
  y0=astr_struc$CRPIX[2]
  xscale=astr_struc$CD[1,1]
  yscale=astr_struc$CD[2,2]
  if (diagnostic) { message(paste("ASTR:",ra0,dec0,x0,y0,xscale,yscale)) }
  pix<-radec2xy(ra,dec,ra0,dec0,x0,y0,xscale,yscale)
  return(pix)
}

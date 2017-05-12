ad.to.xy <-
function(ra,dec,astr.struc,diagnostic=FALSE){
#Details {{{
# Converts RA/DEC (degrees) to x/y (pixels) position using the TAN Gnomonic projection system
# point of interest RA & DEC, anchor point RA & DEC, anchor point x & y, x & y scale (degrees per pixel)
# FROM: http://mathworld.wolfram.com/GnomonicProjection.html
#}}}
  if (missing(dec) & length(dim(ra))==2) { 
    if (ncol(ra)!=2) { stop(paste0('ra must be of dim 2, not ',ncol(ra))) }
    dec<-ra[,2]
    ra<-ra[,1]
  }
  #Get Astrometry Values {{{
  ra0=astr.struc$CRVAL[1]
  dec0=astr.struc$CRVAL[2]
  x0=astr.struc$CRPIX[1]
  y0=astr.struc$CRPIX[2]
  xscale=astr.struc$CD[1,1]
  yscale=astr.struc$CD[2,2]
  #}}}
  #Diagnostic {{{
  if (diagnostic) { message(paste("ASTR:",ra0,dec0,x0,y0,xscale,yscale)) }
  #}}}
  #Get Pixel Values & Return {{{
  pix<-radec.to.xy(ra,dec,ra0,dec0,x0,y0,xscale,yscale)
  return=pix
  #}}}
}

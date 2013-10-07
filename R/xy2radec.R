xy2radec <-
function(x,y,ra0=0,dec0=0,x0=0,y0=0,xscale=1,yscale=1) {
  # Converts x/y (pixels) to RA/DEC (degrees) position using the TAN Gnomonic projection system
  # point of interest x & y, anchor point RA & DEC, anchor point x & y, x & y scale (degrees per pixel)
  # FROM: http://mathworld.wolfram.com/GnomonicProjection.html
  ra0=ra0*(pi/180)
  dec0=dec0*(pi/180)
  xscale=xscale*(pi/180)
  yscale=yscale*(pi/180)
  x = (x-x0)*tan(xscale)
  y = (y-y0)*tan(yscale)
  xra = function(ra0,dec0,x,y){
      ra0 + atan2(x*sin(atan(sqrt(x^2+y^2))),sqrt(x^2+y^2)*cos(dec0)*cos(atan(sqrt(x^2+y^2))) - y*sin(dec0)*sin(atan(sqrt(x^2+y^2))))
  }
  ydec = function(dec0,x,y){
      asin(cos(atan(sqrt(x^2+y^2)))*sin(dec0) + (y*sin(atan(sqrt(x^2+y^2)))*cos(dec0) / sqrt(x^2+y^2))) 
  }
  RA = xra(ra0,dec0,x,y)*180/pi
  DEC = ydec(dec0,x,y)*180/pi
  DEC[which(is.nan(DEC))] = dec0*180/pi
  return(cbind(RA,DEC))
}

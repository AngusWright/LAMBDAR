radec2xy <-
function(ra,dec,ra0=0,dec0=0,x0=0,y0=0,xscale=1,yscale=1){
# Converts RA/DEC (degrees) to x/y (pixels) position using the TAN Gnomonic projection system
# point of interest RA & DEC, anchor point RA & DEC, anchor point x & y, x & y scale (degrees per pixel)
# FROM: http://mathworld.wolfram.com/GnomonicProjection.html
  ra0=ra0*(pi/180)
  dec0=dec0*(pi/180)
  ra=ra*(pi/180)
  dec=dec*(pi/180)
  xscale=xscale*(pi/180)
  yscale=yscale*(pi/180)
  xx = function(ra0,dec0,ra,dec){
          (cos(dec)*sin(ra-ra0))/(sin(dec0)*sin(dec)+(cos(dec0)*cos(dec)*cos(ra-ra0)))
  }
  yy = function(ra0,dec0,ra,dec){
          ((cos(dec0)*sin(dec))-(sin(dec0)*cos(dec)*cos(ra-ra0)))/(sin(dec0)*sin(dec)+(cos(dec0)*cos(dec)*cos(ra-ra0)))
  }
  X = (xx(ra0,dec0,ra,dec)/tan(xscale)) + x0
  Y = (yy(ra0,dec0,ra,dec)/tan(yscale)) + y0
  return=cbind(X,Y)
}

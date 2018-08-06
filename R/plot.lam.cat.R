plot.lam.cat<-function(par.file='Lambdar_default.par',quiet=FALSE) { 

  environment(open.catalogue)<-environment()
  environment(read.images)<-environment()
  environment(read.par.file)<-environment()
  for (nam in ls.deb("package:LAMBDAR",simple=TRUE)) { 
    if (nam%in%ls(envir=environment())) { 
      debug(get(nam,envir=environment()))
    }
  }
  #Read the parameter file
  param.env<-new.env(parent=environment())
  image.env<-new.env(parent=environment())
  param.warnings<-read.par.file(par.file,0,quiet,env=param.env)
  parameter.list<-ls(envir=param.env)
  get.nth.var(parameter.list,n=1,inenv=param.env,outenv=environment(),lim=1)
  #read the catalogue 
  open.catalogue(outenv=environment(),save.table=FALSE)
  if (!filt.contam) { 
    contams<-rep(FALSE,length(cat.ra))
  }
  #read the image
  image.read<-try(read.images(env=NULL,quiet,FALSE,outenv=image.env))
  
  if (ang.offset) { theta.offset<-90-cat.theta } else {theta.offset<-cat.theta}
  #Correct for any reversal of fits image
  if (image.env$astr.struc$CD[1,1]*image.env$astr.struc$CD[2,2]<0){ theta.offset<-theta.offset*-1 }

  #Plot the Input and Simulation Images
  layout(cbind(1,2))

  if (class(image.read)=='try-error') { 
    #No images, just plot the apertures
    #Plot Input & Aperture Ellipses
    magplot(cat.ra,cat.dec,pch=3,cex=0.5,col=ifelse(contams,'red','green3'),
            xlab='RA (deg)',ylab="DEC (deg)")
    invisible(Vectorize(function(x,y,a,b,pa,col) {
                          lines(ellipse(x,y,a,b,pa=pa),col=col,lty=2)
                          points(x,y,pch=3,cex=0.5,col=col)
                        })(cat.ra,cat.dec,cat.a/3600,cat.b/3600,theta.offset-90,ifelse(contams,'red','green3')))
    #Plot Zoomed Sections of the Image:
    #-> Input & Aperture Ellipses
    magplot(cat.ra,cat.dec,pch=3,cex=0.5,col=ifelse(contams,'red','green3'),
            xlab='RA (deg)',ylab="DEC (deg)",xlim=median(cat.ra)+c(-50,+50)/3600,
            ylim=median(cat.dec)+c(-50,+50)/3600)
    invisible(Vectorize(function(x,y,a,b,pa,col) {
                          lines(ellipse(x,y,a,b,pa=pa),col=col,lty=2)
                          points(x,y,pch=3,cex=0.5,col=col)
                        })(cat.ra,cat.dec,cat.a/3600,cat.b/3600,theta.offset-90,ifelse(contams,'red','green3')))
  } else { 
    #Get Object Positions
    xy<-ad.to.xy(cat.ra,cat.dec,image.env$astr.struc)
    aspp<-abs(image.env$astr.struc$CD[1,1])*3600
    #Plot Input & Aperture Ellipses
    magimage(image.env$im)
    #image(x=1:length(image.env$im[,1]),y=1:length(image.env$im[1,]),
    #      image.env$im,col=grey.colors(1E3),zlim=quantile(image.env$im,c(0,0.999)),
    #      asp=1,useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)")
    invisible(Vectorize(function(x,y,a,b,pa,col) {
                          lines(ellipse(x,y,a,b,pa=pa),col=col,lty=2)
                          points(x,y,pch=3,cex=0.5,col=col)
                        })(xy[,1]-1.5,xy[,2]-1.5,cat.a/aspp,cat.b/aspp,theta.offset-90,ifelse(contams,'red','green3')))
    #points(xy[,1],xy[,2],pch=3,cex=0.5,col=ifelse(contams,'red','green3'))
    #Plot Zoomed Sections of the Image:
    #-> Input & Aperture Ellipses
    magimage(image.env$im,xlim=xy[which(!contams)[1],1]+c(-100,+100),
                ylim=xy[which(!contams)[1],2]+c(-100,+100))
    #image(x=1:length(image.env$im[,1]),y=1:length(image.env$im[1,]),
    #      image.env$im,col=grey.colors(1E3),zlim=quantile(image.env$im,c(0,0.999)),
    #      asp=1,useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)",
    #      xlim=image.env$astr.struc$NAXIS[1]/2+c(-50,+50),ylim=image.env$astr.struc$NAXIS[2]/2+c(-50,+50))
    invisible(Vectorize(function(x,y,a,b,pa,col) {
                          lines(ellipse(x,y,a,b,pa=pa),col=col,lty=2)
                          points(x,y,pch=3,cex=0.5,col=col)
                        })(xy[,1]-1.5,xy[,2]-1.5,cat.a/aspp,cat.b/aspp,theta.offset-90,ifelse(contams,'red','green3')))
    #points(xy[,1],xy[,2],pch=3,cex=0.5,col=ifelse(contams,'red','green3'))
  } 

} 

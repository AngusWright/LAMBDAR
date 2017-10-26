#
#
# Function to truncate the upturn (from confusion) 
# in an estimate of the PSF 
#
#

truncate.upturn<-function(psf,centre,poly.degree=8,tolerance=0.01,min.rad=3,flexible=TRUE,cutdown=FALSE,plot=FALSE) { 

  #if no centre specified, use the maxima
  if (missing(centre)) { 
    centre<-which(psf==max(psf,na.rm=TRUE),arr.ind=TRUE)
  }
  #Produce CoG 
  cog<-get.cog(psf,centre=centre,poly.degree=poly.degree,flexible=flexible)
  centre<-cog$cen
  if (plot) { 
    magplot(cog$avg$x,cog$avg$y,type='s',xlab='radius (pix)',ylab='Enclosed Integral (pix units)',ylim=range(cog$avg$y[which(cog$avg$y>0)]),log='y') 
    points(cog$all,pch='.')
    lines(cog$avg$x,cog$avg$fitted,type='s',col='grey') 
    legend('bottomright',inset=0.1,bty='n',legend=c('Raw','Smoothed'),lty=1,col=c('black','grey')) 
    magplot(cog$avg$x,cog$avg$slope,type='s',xlab='radius (pix)',ylab='First Derivative') 
    magplot(cog$avg$x,cog$avg$concav,type='s',xlab='radius (pix)',ylab='Second Derivative') 
  } 
  #Calculate upturn point
  upturn.i=which(abs(cog$avg$concav) <= diff(range(cog$avg$concav))*tolerance)
  upturn=cog$avg$x[upturn.i]
  if (length(upturn)==0) { 
    cat("No upturn found!\n") 
  } else { 
    if (plot) { points(upturn,cog$avg$concav[upturn.i],col='red') }
    upturn=density(upturn,weight=rep(1,length(upturn)),bw=diff(range(cog$avg$x))*tolerance/sqrt(12),kern='rect',from=0,n=1e3)
    upturn.i=which(upturn$y>1 & upturn$x > min.rad)
    start=upturn$x[upturn.i[1]]
    end=upturn$x[upturn.i[which((seq(1:length(upturn.i))-1+upturn.i[1])!=upturn.i)[1]-1]]
    upturn=(start+end)/2
    if (plot) { 
      abline(v=upturn,col='blue',lwd=1.5,lty=2) 
      abline(v=start,col='blue',lwd=1.5,lty=3) 
      abline(v=end,col='blue',lwd=1.5,lty=3) 
    } 
    #Truncate pixels with radius > upturn
    x = seq(1,dim(psf)[1])
    y = seq(1,dim(psf)[2])
    xy = expand.grid(x,y)
    r=sqrt((xy[,1]-centre[1])^2+(xy[,2]-centre[2])^2)
    xy<-xy[which(r>=upturn),]
    if (plot) { 
      magimage(psf)
      lines(ellipse(xcen=centre[1]-0.5,ycen=centre[2]-0.5,a=upturn,e=0,pa=0),col='lightblue',lty=2,lwd=1.5)
      lines(ellipse(xcen=centre[1]-0.5,ycen=centre[2]-0.5,a=start,e=0,pa=0),col='lightblue',lty=3,lwd=1.5)
      lines(ellipse(xcen=centre[1]-0.5,ycen=centre[2]-0.5,a=end,e=0,pa=0),col='lightblue',lty=3,lwd=1.5)
    }
    psf[cbind(xy[,1],xy[,2])]<-0
    if (plot) { 
      #Final plots 
      magplot(cog$avg$x,cog$avg$y,type='s',xlab='radius (pix)',ylab='Enclosed Integral (pix units)',col='grey',ylim=range(cog$avg$y[which(cog$avg$y>0)]),log='y') 
      points(cog$all,pch='.',col='grey')
      abline(v=upturn,col='blue',lwd=1.5,lty=2) 
      cog<-get.cog(psf,centre=centre,flexible=FALSE)
      points(cog$all,pch='.') 
      lines(cog$avg$x,cog$avg$y,type='s') 
      magimage(psf)
    } 
    if (cutdown) { 
      xy = expand.grid(x,y)
      xy = xy[which(r<=upturn),]
      xin= range(xy[,1])
      yin= range(xy[,2])
      psf<-psf[xin[1]:xin[2],yin[1]:yin[2]]
    }
  }
  return=psf
}


para.read.fits<-
  function(filename,nsubim=2,outputfname="crop"){
 
  #Read Header to determine number of pixels
  bitpix<-abs(as.numeric(read.fitskey("BITPIX",filename)))
  #Read image axis lengths
  naxisx<-as.numeric(read.fitskey("NAXIS1",filename))
  naxisy<-as.numeric(read.fitskey("NAXIS2",filename))
  rat<-naxisx/naxisy
  #Determine number of pixels wanted in sub-images
  npix<-(naxisx*naxisy)/nsubim
  if (npix-floor(npix) != 0) {
    for(i in 1:nsubim) { if (i%%2!=0) { nspix<-c(nspix,floor(npix)) } else { nspix<-c(nspix, ceiling(npix)) } }
  }
  if (sum(nspix)!=(naxisx*naxisy)) { nspix[nsubim]<-nspix[nsubim]+((naxisx*naxisy)-sum(nspix)) }
  subimsize<-npix*bitpix
  #Determine number of sub-images we have
  nysub<-round(nsubim/rat, digits=0)
  nxsub<-round(nsubim*rat, digits=0)
  if (nysub*nxsub != nsubim) { nxsub<-nsubim-(nysub*nxsub)/nysub }
  #Determine sub-image axis lengths
  nypix<-sqrt(nspix/rat)
  if (max(nypix-floor(nypix)) != 0) {
    for(i in 1:nysub) { if (i%%2!=0) { nsypix<-c(nsypix,floor(nypix)) } else { nsypix<-c(nsypix, ceiling(nypix)) } }
  }
  if (sum(nsypix)!=sum(nypix)) { nsypix[nsubim]<-nsypix[nsubim]+(sum(nypix)-sum(nsypix)) }
  nxpix<-sqrt(nspix*rat)
  if (max(nxpix-floor(nxpix)) != 0) {
    for(i in 1:nxsub) { if (i%%2!=0) { nsxpix<-c(nsxpix,floor(nxpix)) } else { nsxpix<-c(nsxpix, ceiling(nxpix)) } }
  }
  if (sum(nsxpix)!=sum(nxpix)) { nsxpix[nsubim]<-nsxpix[nsubim]+(sum(nxpix)-sum(nsxpix)) }
  #make pixel indicies cumulative
  xsum<-1
  xpixlim<-NULL
  for (i in 1:nxsub) { 
    xpixlim<-rbind(xpixlim, c(xsum,xsum+nsxpix[i]))
    xsum<-xsum+nsxpix[i] 
  }
  ysum<-1
  ypixlim<-NULL
  for (i in 1:nysub) { 
    ypixlim<-rbind(ypixlim, c(ysum,ysum+nsypix[i]))
    ysum<-ysum+nsypix[i] 
  }
  message(paste("There are",nxsub*nysub,"sub-images to read"))
  message(paste("SubImages are roughly",subimsize,"bits in size"))
  message(paste("ImAxis Lengths are",naxisx,"and",naxisy))
  message(paste("SubAxis Lengths are",nxpix,"and",nypix))
  message(paste("SubImages Indicies are",nxsub,"and",nysub))
  #Remove any previous subimages
  system(paste("rm -f ",outputfname,"*", sep=""))
  #Get the subimages
  im<-foreach(j=1:nysub, yaxl=ypixlim[,1], yaxu=ypixlim[,2], .combine='rbind')%:%
      foreach(i=1:nxsub, xaxl=xpixlim[,1], xaxu=xpixlim[,2], .combine='cbind')%dopar%{
    #Make cropped image
    cutcommandim <- paste("fitscopy ",filename,"[",format(xaxl,scientific=FALSE),":",format(xaxu,scientific=FALSE),",",
						   format(yaxl,scientific=FALSE),":",format(yaxu,scientific=FALSE),"] ",
						   outputfname,"_",i,"_",j,".fits",sep="")
    system(cutcommandim)
    #Read cropped image 
    sub<-read.fits(paste(outputfname,"_",i,"_",j,".fits",sep=""),hdu=0,comment=FALSE)
    sub$dat[[1]]
  }
  return(im)
}

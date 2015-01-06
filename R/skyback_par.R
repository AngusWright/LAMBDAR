skyback.par<-function(x_p,y_p,cutlo=0,cuthi=100,im_mask,imm_mask,remmask=TRUE,radweight=1,clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,probcut=3, mpiopts=""){
  if(length(x_p) != length(y_p)){stop('x_pix and y_pix lengths do not much!')}
  if(length(x_p) != length(cutlo)){stop('x_pix and cutlo lengths do not much!')}
  if(length(x_p) != length(cuthi)){stop('x_pix and cuthi lengths do not much!')}
  if(length(x_p) != length(im_mask)){stop('x_pix and im_mask lengths do not much!')}
  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(probcut)

  output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x_p, pixlocy=y_p, origim=im_mask, maskim=imm_mask, .options.mpi=mpiopts, .combine='rbind') %dopar% {
  #browser()
  #for (kk in 346:length(cutlovec)) {
  #  cutlo=cutlovec[kk]
  #  cuthi=cuthivec[kk]
  #  pixlocx=x_p[kk]
  #  pixlocy=y_p[kk]
  #  origim=im_mask[[kk]]
  #  maskim=imm_mask[[kk]]

    #Find extreme pixels to cut out
    xlocs<-pixlocx+(-cuthi:cuthi)
    ylocs<-pixlocy+(-cuthi:cuthi)
    #Find object location on new pixel grid
    xcen<-cuthi+1-length(which(xlocs<0))
    ycen<-cuthi+1-length(which(ylocs<0))
    #Select only pixels which are inside the image bounds
    xsel<-which(xlocs>0 & xlocs<=length(origim[,1]))
    ysel<-which(ylocs>0 & ylocs<=length(origim[1,]))
    #Trim to above
    xlocs<-xlocs[xsel]
    ylocs<-ylocs[ysel]
    #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
    if(remmask){
      tempim<-origim*maskim
    }else{
      tempim<-origim
    }
    #All ref pixels for new image
    tempref<-as.matrix(expand.grid(1:length(xsel),1:length(ysel)))
    #Corresponding radii for new pixels from the object of interest
    temprad<-sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
    #Keep only pixels inside the radius bounds given by cutlo and cuthi
    keep<-temprad>cutlo & temprad<cuthi
    #Trim
    tempref<-tempref[keep,]
    tempval<-tempim[tempref]
    temprad<-temprad[keep]
    #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
    if(remmask){
      temprad<-temprad[tempval!=0]
      tempval<-tempval[tempval!=0]
    }
    #Do iterative <probcut>-sigma pixel clipping
    if(clipiters>0){
      for(j in 1:clipiters){
        vallims<-2*median(tempval)-quantile(tempval,probcut, na.rm=TRUE)
        temprad<-temprad[tempval<vallims]
        tempval<-tempval[tempval<vallims]
      }
    }
    #Find the running medians for the data
    tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T)
    tempylims<-tempmedian$ysd
    tempy<-tempmedian$y
    tempx<-tempmedian$x
    #Remove bins with no skypixels present
    if (any(is.na(tempylims))) {
      tempy  <-tempy[which(!is.na(tempylims[,1]))]
      tempx  <-tempx[which(!is.na(tempylims[,1]))]
      temprad<-temprad[which(!is.na(tempylims[,1]))]
      tempval<-tempval[which(!is.na(tempylims[,1]))]
      tempref<-tempref[which(!is.na(tempylims[,1]))]
      tempylims<-matrix(tempylims[which(!is.na(tempylims),arr.ind=TRUE)], ncol=2)
    }
    if (length(tempy)!=0) {
      #Calculate worst case sky error- the sd of the medians calculated
      skyerr<-sd(tempy)
      #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
      weights<-1/((tempx^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
      #Generate Sky RMS
      skyRMS<-as.numeric((quantile(tempval,0.5, na.rm=TRUE)-quantile(tempval,pnorm(-1),na.rm=TRUE)))
      #Determine Pearsons Test for Normality p-value for sky
      skyRMSpval<-pearson.test(tempval)$p.value
      #Find the weighted mean of the medians
      sky<-sum(tempy*weights)/(sum(weights))
      #Now we iterate until no running medians are outside the 1-sigma bound of the sky
      nloop<-0
      while(any(!(tempylims[,1]<=sky & tempylims[,2]>=sky)) & all(!(tempylims[,1]<=sky & tempylims[,2]>=sky))==FALSE){
        nloop<-nloop+1
        tempy<-tempy[tempylims[,1]<=sky & tempylims[,2]>=sky]
        weights<-weights[tempylims[,1]<=sky & tempylims[,2]>=sky]
        tempylims<-rbind(tempylims[tempylims[,1]<=sky & tempylims[,2]>=sky,])
        sky<-sum(tempy*weights)/(sum(weights))
        if (nloop==1E8) {
          warning("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
          message("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
          break
        }
      }
      #Find the number of running medians that agree with the final sky within error bounds (max=10)
      Nnearsky<-length(which(tempylims[,1]<=sky & tempylims[,2]>=sky))
      #Organise data into data.frame for foreach
      return=data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval)
    } else {
      #No Sky estimate available - return NAs
      return=data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA)
    }
  }
  #Output the foreach data
  return=output
}

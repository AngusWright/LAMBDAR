skyback<-function(ra,dec,cutlo=0,cuthi=100,origim,astrom,maskim,remmask=TRUE,radweight=1,clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,probcut=3){
  if(length(ra) != length(dec)){stop('ra and dec lengths do not much!')}
  if(length(ra) != length(cutlo)){stop('ra and cutlo lengths do not much!')}
  if(length(ra) != length(cuthi)){stop('ra and cuthi lengths do not much!')}
  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  ravec<-ra;decvec<-dec;cutlovec<-round(cutlo);cuthivec<-round(cuthi)
  pixlocvec<-round(ad2xy(ra,dec,astrom))
  count<-length(ra)
  probcut<-1-pnorm(probcut)

  output<-NULL
  estimated<-FALSE
  for (i in 1:count){
    ra<-ravec[i];dec<-decvec[i];cutlo<-cutlovec[i];cuthi<-cuthivec[i];pixloc<-pixlocvec[i,]
    #Find extreme pixels to cut out
    xlocs<-pixloc[1]+(-cuthi:cuthi)
    ylocs<-pixloc[2]+(-cuthi:cuthi)
    #Find object location on new pixel grid
    xcen<-cuthi+1-length(which(xlocs<0))
    ycen<-cuthi+1-length(which(ylocs<0))
    #Select only pixels which are inside the image bounds
    xsel<-which(xlocs>0 & xlocs<=astrom$NAXIS[1])
    ysel<-which(ylocs>0 & ylocs<=astrom$NAXIS[2])
    #Trim to above
    xlocs<-xlocs[xsel]
    ylocs<-ylocs[ysel]
    #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
    if(remmask){
      tempim<-origim$dat[[1]][xlocs,ylocs]*maskim$dat[[1]][xlocs,ylocs]
    }else{
      tempim<-origim$dat[[1]][xlocs,ylocs]
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
      for(i in 1:clipiters){
        vallims<-2*median(tempval)-quantile(tempval,probcut)
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
      tempy  <-tempy[which(!is.na(tempylims))]
      tempx  <-tempx[which(!is.na(tempylims))]
      temprad<-temprad[which(!is.na(tempylims))]
      tempval<-tempval[which(!is.na(tempylims))]
      tempref<-tempref[which(!is.na(tempylims))]
    }
    if (length(tempy)!=0) {
      #Calculate worst case sky error- the sd of the medians calculated
      skyerr<-sd(tempy)
      #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
      weights<-1/((tempx^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
      #Generate Sky RMS
      skyRMS<-as.numeric((quantile(tempval,0.5)-quantile(tempval,pnorm(-1))))
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
      output<-rbind(output,data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval))
    } else {
      #No Sky estimate available - return NAs
      output<-rbind(output,data.frame(sky=0,skyerr=0,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA))
    }
  }
  gc()
  #Output the foreach data
  return=output
}

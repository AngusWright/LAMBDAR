sky.estimate<-function(x.pix,y.pix,cutlo=0,cuthi=100,data.stamp,mask.stamp=NULL,remmask=TRUE,radweight=1,clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,probcut=3, mpi.opts=""){
  cutup=TRUE
  if(length(x.pix) != length(y.pix)){stop('x.pixix and y.pixix lengths do not match!')}
  if(length(x.pix) != length(cutlo)){stop('x.pixix and cutlo lengths do not match!')}
  if(length(x.pix) != length(cuthi)){stop('x.pixix and cuthi lengths do not match!')}
  if(length(x.pix) != length(data.stamp)){
    if (is.matrix(data.stamp)) {
      cutup=FALSE
    } else {
      stop('x.pixix and data.stamp lengths do not match!')
    }
  }
  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(probcut)

  if (cutup) {
    output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x.pix, pixlocy=y.pix, origim=data.stamp, maskim=mask.stamp, .options.mpi=mpi.opts, .combine='rbind') %dopar% {

      for (run in 1:2) {
        #Find extreme pixels to cut out
        xlocs<-floor(pixlocx+(-cuthi:cuthi))
        ylocs<-floor(pixlocy+(-cuthi:cuthi))
        #Find object location on new pixel grid
        xcen<-cuthi+1-length(which(xlocs<0))
        ycen<-cuthi+1-length(which(ylocs<0))
        #Select only pixels which are inside the image bounds
        xsel<-which(xlocs>0 & xlocs<=length(origim[,1]))
        ysel<-which(ylocs>0 & ylocs<=length(origim[1,]))
        #Trim to above
        xlocs<-xlocs[xsel]
        ylocs<-ylocs[ysel]
        #All ref pixels for new image
        tempref<-as.matrix(expand.grid(1:length(xsel),1:length(ysel)))
        #Corresponding radii for new pixels from the object of interest
        temprad<-sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
        #Keep only pixels inside the radius bounds given by cutlo and cuthi
        keep<-temprad>cutlo & temprad<cuthi
        #Trim
        tempref<-tempref[keep,]
        if(remmask){
          tempval<-origim[tempref]*maskim[tempref]
        }else{
          tempval<-origim[tempref]
        }
        temprad<-temprad[keep]
        #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
        if(remmask){
          temprad<-temprad[tempval!=0]
          tempval<-tempval[tempval!=0]
        }
        #Do iterative <probcut>-sigma pixel clipping
        if(clipiters>0){
          for(j in 1:clipiters){
            if (run==1) {
              vallims<-2*median(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            } else {
              vallims<-2*mean(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            }
            temprad<-temprad[tempval<vallims]
            tempval<-tempval[tempval<vallims]
          }
        }
        #Find the running medians(run1) or means(run2) for the data
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE,type='mean') }
        tempylims<-tempmedian$ysd
        tempy<-tempmedian$y
        tempx<-tempmedian$x
        #Remove bins with 1 or less skypixels present
        if (any(is.na(tempylims)|(tempylims[,2]==tempylims[,1]))) {
          ind<-which((!is.na(tempylims[,1]))&(!is.na(tempylims[,2]))&(tempylims[,2]!=tempylims[,1]))
          tempy  <-tempy[ind]
          tempx  <-tempx[ind]
          tempylims<-rbind(tempylims[ind,])
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
          if (!is.na(skyerr)) {
            #Now we iterate until no running medians are outside the 1-sigma bound of the sky
            nloop<-0
            while(any(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr))) & all(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))==FALSE){
              nloop<-nloop+1
              ref<-which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr))
              tempy<-tempy[ref]
              weights<-weights[ref]
              tempylims<-rbind(tempylims[ref,])
              sky<-sum(tempy*weights)/(sum(weights))
              if (nloop>15) {
                warning("Failure to converge to within the 1-sigma bound of the sky. Breaking.")
                message("Failure to converge to within the 1-sigma bound of the sky. Breaking.")
                break
              }
            }
            #Find the number of running medians that agree with the final sky within error bounds (max=10)
            Nnearsky<-length(which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))
          } else {
            #Only 1 bin is left
            Nnearsky<- 1
          }
          #Organise data into data.frame for foreach
          if (run==1) { medianStat<-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
          if (run==2) { meanStat  <-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
        } else {
          #No Sky estimate available - return NAs
          if (run==1) { medianStat<-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
          if (run==2) { meanStat  <-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
        }
      }
      return=data.frame(sky=medianStat$sky,skyerr=medianStat$skyerr,Nnearsky=medianStat$Nnearsky,skyRMS=medianStat$skyRMS,skyRMSpval=medianStat$skyRMSpval,
                        sky.mean=meanStat$sky,skyerr.mean=meanStat$skyerr,Nnearsky.mean=meanStat$Nnearsky,skyRMS.mean=meanStat$skyRMS,skyRMSpval.mean=meanStat$skyRMSpval)
    }
  } else {
    output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x.pix, pixlocy=y.pix, .export=c('data.stamp','mask.stamp'), .options.mpi=mpi.opts, .combine='rbind') %dopar% {

      for (run in 1:2) {
        #Find extreme pixels to cut out
        xlocs<-floor(pixlocx+(-cuthi:cuthi))
        ylocs<-floor(pixlocy+(-cuthi:cuthi))
        #Find object location on new pixel grid
        xcen<-cuthi+1-length(which(xlocs<0))
        ycen<-cuthi+1-length(which(ylocs<0))
        #Select only pixels which are inside the image bounds
        xsel<-which(xlocs>0 & xlocs<=length(data.stamp[,1]))
        ysel<-which(ylocs>0 & ylocs<=length(data.stamp[1,]))
        #Trim to above
        xlocs<-xlocs[xsel]
        ylocs<-ylocs[ysel]
        #All ref pixels for new image
        tempref<-as.matrix(expand.grid(xlocs,ylocs))
        #Corresponding radii for new pixels from the object of interest
        temprad<-sqrt((tempref[,1]-pixlocx)^2+(tempref[,2]-pixlocy)^2)
        #Keep only pixels inside the radius bounds given by cutlo and cuthi
        keep<-temprad>cutlo & temprad<cuthi
        #Trim
        tempref<-tempref[keep,]
        #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
        if(remmask){
          tempval<-data.stamp[tempref]*mask.stamp[tempref]
        }else{
          tempval<-data.stamp[tempref]
        }
        temprad<-temprad[keep]
        #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
        if(remmask){
          temprad<-temprad[tempval!=0]
          tempval<-tempval[tempval!=0]
        }
        #Do iterative <probcut>-sigma pixel clipping
        if(clipiters>0){
          for(j in 1:clipiters){
            if (run==1) {
              vallims<-2*median(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            } else {
              vallims<-2*mean(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            }
            temprad<-temprad[tempval<vallims]
            tempval<-tempval[tempval<vallims]
          }
        }
        #Find the running medians(run1) or means(run2) for the data
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE,type='mean') }
        tempylims<-tempmedian$ysd
        tempy<-tempmedian$y
        tempx<-tempmedian$x
        #Remove bins with 1 or less skypixels present
        if (any(is.na(tempylims)|(tempylims[,2]==tempylims[,1]))) {
          ind<-which((!is.na(tempylims[,1]))&(!is.na(tempylims[,2]))&(tempylims[,2]!=tempylims[,1]))
          tempy  <-tempy[ind]
          tempx  <-tempx[ind]
          tempylims<-rbind(tempylims[ind,])
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
          if (!is.na(skyerr)) {
            #Now we iterate until no running medians are outside the 1-sigma bound of the sky
            nloop<-0
            while(any(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr))) & all(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))==FALSE){
              nloop<-nloop+1
              ref<-which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr))
              tempy<-tempy[ref]
              weights<-weights[ref]
              tempylims<-rbind(tempylims[ref,])
              sky<-sum(tempy*weights)/(sum(weights))
              if (nloop>15) {
                warning("Failure to converge to within the 1-sigma bound of the sky. Breaking.")
                message("Failure to converge to within the 1-sigma bound of the sky. Breaking.")
                break
              }
            }
            #Find the number of running medians that agree with the final sky within error bounds (max=10)
            Nnearsky<-length(which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))
          } else {
            #Only 1 bin is left
            Nnearsky<- 1
          }
          #Organise data into data.frame for foreach
          if (run==1) { medianStat<-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
          if (run==2) { meanStat  <-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
        } else {
          #No Sky estimate available - return NAs
          if (run==1) { medianStat<-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
          if (run==2) { meanStat  <-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
        }
      }
      return=data.frame(sky=medianStat$sky,skyerr=medianStat$skyerr,Nnearsky=medianStat$Nnearsky,skyRMS=medianStat$skyRMS,skyRMSpval=medianStat$skyRMSpval,
                        sky.mean=meanStat$sky,skyerr.mean=meanStat$skyerr,Nnearsky.mean=meanStat$Nnearsky,skyRMS.mean=meanStat$skyRMS,skyRMSpval.mean=meanStat$skyRMSpval)
    }
  }
  #Output the foreach data
  return=output
}

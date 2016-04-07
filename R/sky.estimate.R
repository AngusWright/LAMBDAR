#
#
#
sky.estimate<-function(data.stamp,mask.stamp=NULL,rem.mask=FALSE,data.stamp.lims=NULL,cutlo=0,cuthi=100,radweight=1,clipiters=5,PSFFWHMinPIX=0/0.339,hardlo=3,hardhi=10,sigma.cut=3,cat.x=NULL,cat.y=NULL,mpi.opts=NULL){
  #Check Data stamp Limits
  if (is.null(data.stamp.lims)){
    #Data stamps must be the same dimension as the image!
    warning('data.stamp.lims is not provided, so we assume the limits are the stamp edges')
    if (is.matrix(data.stamp)) {
      data.stamp.lims<-matrix(c(1,length(data.stamp[,1]),1,length(data.stamp[1,])),nrow=max(1,length(cat.x)),ncol=4,byrow=T)
    } else {
      data.stamp.lims<-rbind(foreach(ap=data.stamp,.combine='rbind')%dopar%{return=c(1,length(ap[,1]),1,length(ap[1,]))})
    }
  } else if (length(data.stamp) != length(data.stamp.lims[,1])){
    stop('data.stamp.lims is not the same length as data.stamp!')
  }
  if (rem.mask) {
    if (is.null(mask.stamp)) {
      stop("No mask supplied with rem.mask=TRUE!")
    }
    if (any(dim(mask.stamp)!=dim(data.stamp))) {
      stop("Dimensions of mask.stamp and data.stamp are not the same")
    }
  }
  cutup<-TRUE
  if (is.matrix(data.stamp)) {
    cutup<-FALSE
  }

  if (is.null(cat.x)) {
    x.pix<-rowMeans(rbind(data.stamp.lims[,cbind(1,2)]))
  } else {
    x.pix<-floor(cat.x)
  }
  if (is.null(cat.y)) {
    y.pix<-rowMeans(rbind(data.stamp.lims[,cbind(3,4)]))
  } else {
    y.pix<-floor(cat.y)
  }
  if (length(cutlo)!=length(x.pix)) {
    cutlo<-rep(cutlo,length(x.pix))
  }
  if (length(cuthi)!=length(x.pix)) {
    cutlo<-rep(cuthi,length(x.pix))
  }

  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(sigma.cut)

  if (cutup) {
    output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x.pix, pixlocy=y.pix, origim=data.stamp, maskim=mask.stamp, sxl=data.stamp.lims[,1],syl=data.stamp.lims[,3],.options.mpi=mpi.opts, .combine='rbind') %dopar% {

      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
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
        if(rem.mask){
          maskim[which(maskim==0)]<-NA
          tempval<-origim[tempref]*maskim[tempref]
        }else{
          tempval<-origim[tempref]
        }
        temprad<-temprad[keep]
        #Do iterative <probcut> pixel clipping
        if(clipiters>0){
          for(j in 1:clipiters){
            if (run==1) {
              vallims<-2*median(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            } else {
              vallims<-2*mean(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            }
            temprad<-temprad[tempval<=vallims]
            tempval<-tempval[tempval<=vallims]
          }
        }
        #Find the running medians(run1) or means(run2) for the data
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE,type='mean') }
        tempylims<-tempmedian$ysd
        tempy<-tempmedian$y
        tempx<-tempmedian$x
        num.in.bin<-tempmedian$Nbins
        #Remove bins with 1 or less skypixels present
        if (any(num.in.bin<=1)) {
          ind<-which(num.in.bin<=1)
          tempy  <-tempy[ind]
          tempx  <-tempx[ind]
          templty<-templty[ind]
          tempylims<-rbind(tempylims[ind,])
        }
        if (length(tempy)!=0) {
          #Calculate worst case sky error- the sd of the medians calculated
          skyerr<-sd(tempy)
          #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
          #Catch divide by 0
          bin.sd<-(tempylims[,2]-tempylims[,1])/2
          bin.sd[which(bin.sd==0)]<-1E-32
          weights<-1/((tempx^radweight)*bin.sd)^2
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
    output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x.pix, pixlocy=y.pix, sxl=data.stamp.lims[,1],syl=data.stamp.lims[,3], .export=c('data.stamp','mask.stamp'), .options.mpi=mpi.opts, .combine='rbind') %dopar% {

      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
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
        if(rem.mask){
          mask.stamp[which(mask.stamp==0)]<-NA
          tempval<-data.stamp[tempref]*mask.stamp[tempref]
        }else{
          tempval<-data.stamp[tempref]
        }
        temprad<-temprad[keep]
        #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
        #if(rem.mask){
        #  temprad<-temprad[tempval!=0]
        #  tempval<-tempval[tempval!=0]
        #}
        #Do iterative <probcut>-sigma pixel clipping
        if(clipiters>0){
          for(j in 1:clipiters){
            if (run==1) {
              vallims<-2*median(tempval,na.rm=TRUE)-quantile(tempval,probcut, na.rm=TRUE)
            } else {
              vallims<-2*mean(tempval,na.rm=TRUE)-quantile(tempval,probcut, na.rm=TRUE)
            }
            temprad<-temprad[tempval<=vallims]
            tempval<-tempval[tempval<=vallims]
          }
        }
        #Find the running medians(run1) or means(run2) for the data
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE,type='mean') }
        tempylims<-tempmedian$ysd
        tempy<-tempmedian$y
        tempx<-tempmedian$x
        num.in.bin<-tempmedian$Nbins
        #Remove bins with 1 or less skypixels present
        if (any(num.in.bin<=1)) {
          ind<-which(num.in.bin<=1)
          tempy  <-tempy[ind]
          tempx  <-tempx[ind]
          templty<-templty[ind]
          tempylims<-rbind(tempylims[ind,])
        }
        if (length(tempy)!=0) {
          #Calculate worst case sky error- the sd of the medians calculated
          skyerr<-sd(tempy)
          #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
          #Catch divide by 0
          bin.sd<-(tempylims[,2]-tempylims[,1])/2
          bin.sd[which(bin.sd==0)]<-1E-32
          weights<-1/((tempx^radweight)*bin.sd)^2
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

#
#
#
fast.sky.estimate<-function(data.stamp,mask.stamp=NULL,rem.mask=TRUE,data.stamp.lims=NULL,cutlo=0,cuthi=100,radweight=1,clipiters=5,PSFFWHMinPIX=0,hardlo=3,hardhi=10,sigma.cut=3,saturate=Inf,cat.x=NULL,cat.y=NULL,mpi.opts=NULL,subset,bins=10,fit.gauss=FALSE){
  #Check Data stamp Limits
  if (is.null(data.stamp.lims)){
    #Data stamps must be the same dimension as the image!
    warning('data.stamp.lims is not provided, so we assume the limits are the stamp edges')
    if (is.matrix(data.stamp)) {
      data.stamp.lims<-matrix(c(1,length(data.stamp[,1]),1,length(data.stamp[1,])),nrow=max(1,length(cat.x)),ncol=4,byrow=T)
    } else {
      data.stamp.lims<-rbind(foreach(ap=data.stamp,.combine='rbind')%dopar%{return=c(1,length(ap[,1]),1,length(ap[1,]))})
    }
  } else if (!is.matrix(data.stamp) & (length(data.stamp) != length(data.stamp.lims[,1]))){
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
    cuthi<-rep(cuthi,length(x.pix))
  }

  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(sigma.cut)

  if (!missing(subset)) { 
    cutlovec<-cutlovec[subset]
    cuthivec<-cuthivec[subset]
    x.pix<-x.pix[subset]
    y.pix<-y.pix[subset]
    data.stamp.lims<-rbind(data.stamp.lims[subset,])
    if (cutup) { 
      data.stamp<-data.stamp[subset]
      mask.stamp<-mask.stamp[subset]
    }
  } else { 
    subset<-1:length(x.pix) 
  }
  if (length(subset)!=0) { 

    if (cutup) {
      output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x.pix, pixlocy=y.pix, origim=data.stamp, maskim=mask.stamp, sxl=data.stamp.lims[,1],syl=data.stamp.lims[,3],.options.mpi=mpi.opts, .combine='rbind') %dopar% {

        pixlocx=pixlocx-(sxl-1)
        pixlocy=pixlocy-(syl-1)
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
          tempval<-origim[tempref][which(maskim[tempref]!=0)]
        }else{
          tempval<-origim[tempref]
        }
        #Remove saturated pixels 
        if(length(tempval)!=0) { 
          tempval<-tempval[which(tempval < saturate)]
        }
        #Do iterative <probcut> pixel clipping
        if(length(tempval)!=0 & clipiters>0){
          for(j in 1:clipiters){
            vallims<-2*median(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            tempval<-tempval[which(tempval<=vallims)]
          }
        }
        #Fit gaussian to the low end of the data
        if (length(tempval)!=0 & fit.gauss) { 
          fit<-fit.gauss2low(tempval,na.rm=T)
        } else { 
          fit<-data.frame(mu=NA,muerr=NA,sd=NA)
        }
        
        rms<-as.numeric((quantile(tempval,0.5, na.rm=TRUE)-quantile(tempval,pnorm(-1),na.rm=TRUE)))
        return=data.frame(skyMu=fit$mu,skyMuErr=fit$muerr,Nnearsky=NA,skySD=fit$sd,skyRMSpval=pearson.test(tempval)$p.value,
                          skyMedian=median(tempval),skyMedianErr=rms/sqrt(length(tempval)),skyRMS=rms)
      }
    } else {
      output<-foreach(cutlo=cutlovec, cuthi=cuthivec, pixlocx=x.pix, pixlocy=y.pix, sxl=data.stamp.lims[,1],syl=data.stamp.lims[,3], .export=c('data.stamp','mask.stamp'), .options.mpi=mpi.opts, .combine='rbind') %dopar% {

        pixlocx=pixlocx-(sxl-1)
        pixlocy=pixlocy-(syl-1)
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
          tempval<-data.stamp[tempref][which(mask.stamp[tempref]!=0)]
        }else{
          tempval<-data.stamp[tempref]
        }
        #Remove saturated pixels 
        if(length(tempval)!=0) { 
          tempval<-tempval[which(tempval < saturate)]
        }
        #Do iterative <probcut> pixel clipping
        if(length(tempval)!=0 & clipiters>0){
          for(j in 1:clipiters){
            vallims<-2*median(tempval)-quantile(tempval,probcut, na.rm=TRUE)
            tempval<-tempval[which(tempval<=vallims)]
          }
        }
        #Fit gaussian to the low end of the data
        if (length(tempval)!=0 & fit.gauss) { 
          fit<-fit.gauss2low(tempval,na.rm=T)
        } else { 
          fit<-data.frame(mu=NA,muerr=NA,sd=NA)
        }
        
        rms<-as.numeric((quantile(tempval,0.5, na.rm=TRUE)-quantile(tempval,pnorm(-1),na.rm=TRUE)))
        return=data.frame(skyMu=fit$mu,skyMuErr=fit$muerr,Nnearsky=NA,skySD=fit$sd,skyRMSpval=pearson.test(tempval)$p.value,
                          skyMedian=median(tempval),skyMedianErr=rms/sqrt(length(tempval)),skyRMS=rms)
      }
    }
  } else { 
    output=data.frame(skyMu=NA,skyMuErr=NA,Nnearsky=NA,skySD=NA,skyRMSpval=NA,
                          skyMedian=NA,skyMedianErr=NA,skyRMS=NA)[-1,]
  }
  #Output the foreach data
  return=output
}

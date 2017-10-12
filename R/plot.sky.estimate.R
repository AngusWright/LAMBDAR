#
#
#
plot.sky.estimate<-function(data.stamp,mask.stamp=NULL,rem.mask=FALSE,data.stamp.lims=NULL,cutlo=0,cuthi=100,radweight=1,clipiters=5,PSFFWHMinPIX=0,hardlo=3,hardhi=10,sigma.cut=3,plot.all=FALSE,path=NULL,toFile=FALSE,cat.id=NULL,cat.x=NULL,cat.y=NULL,res=220,bin.lwd=1,est.lwd=1.5,all.est.lwd=1.5){
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
  if (toFile) {
    if (is.null(cat.id)) { stop("Output to file requires cat.id to be specified (used for filenames)") }
    dir.create(file.path(path,"SkyBackIms"),showWarnings=FALSE)
    path=file.path(path,"SkyBackIms")
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
  if (length(cat.id)!=length(x.pix)) {
    cat.id<-rep(NA,length(x.pix))
  }

  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(sigma.cut)

  if(plot.all) {
    rand<-1:length(x.pix)
  } else {
    rand<-sample(length(x.pix),min(10,length(x.pix)))
  }
  if (cutup) {
    for (r in rand) {
      cutlo=cutlovec[r]
      cuthi=cuthivec[r]
      origim=data.stamp[[r]]
      maskim=mask.stamp[[r]]
      sxl=data.stamp.lims[r,1]
      syl=data.stamp.lims[r,3]
      pixlocx=x.pix[r]
      pixlocy=y.pix[r]
      id=cat.id[r]

      pixloc=c(pixlocx,pixlocy)
      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
      for (run in 2:1) {
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
        #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
        if(rem.mask){
          maskim[which(maskim==0)]<-NA
          tempim<-origim*maskim
        }else{
          tempim<-origim
        }
        #All ref pixels for new image
        tempref<-as.matrix(expand.grid(xlocs,ylocs))
        #Corresponding radii for new pixels from the object of interest
        temprad<-sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
        #Keep only pixels inside the radius bounds given by cutlo and cuthi
        keep<-temprad>cutlo & temprad<cuthi
        #Trim
        tempref<-tempref[keep,]
        tempval<-tempim[tempref]
        temprad<-temprad[keep]
        #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
        #if(rem.mask){
        #  temprad<-temprad[tempval!=0]
        #  tempval<-tempval[tempval!=0]
        #}
        #Do iterative <probcut> pixel clipping
        if(clipiters>0){
          for(j in 1:clipiters){
            if (run==1) {
              vallims<-2*median(tempval,na.rm=TRUE)-quantile(tempval,probcut, na.rm=TRUE)
            } else {
              vallims<-2*mean(tempval,na.rm=TRUE)-quantile(tempval,probcut, na.rm=TRUE)
            }
            temprad<-temprad[which(tempval<=vallims)]
            tempval<-tempval[which(tempval<=vallims)]
          }
        }
        #Find the running medians(run1) or means(run2) for the data
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=TRUE,type='mean') }
        tempylims<-tempmedian$ysd
        tempy<-tempmedian$y
        tempx<-tempmedian$x
        templty<-rep(1:2,ceiling(length(tempx)/2))[1:length(tempx)]
        num.in.bin<-tempmedian$Nbins
        #Remove bins with 1 or less skypixels present
        if (any(num.in.bin<=1)) {
          ind<-which(!num.in.bin<=1)
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
            tempybak<-tempy
            tempxbak<-tempx
            templtybak<-templty
            tempylimsbak<-tempylims
            skybak<-sky
            while(any(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr))) & all(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))==FALSE){
              nloop<-nloop+1
              tempy<-tempy[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              tempx<-tempx[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              templty<-templty[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              weights<-weights[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              tempylims<-rbind(tempylims[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr),])
              sky<-sum(tempy*weights)/(sum(weights))
              if (nloop>15) {
                warning("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
                message("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
                break
              }
            }
            #Find the number of running medians that agree with the final sky within error bounds (max=10)
            Nnearsky<-length(which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))
          } else {
            Nnearsky<-1
          }
          if (run==1) { medianStat<-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
          if (run==2) { meanStat  <-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
        } else {
          #No Sky estimate available - return NAs
          sky=NA
          skyerr=NA
          Nnearsky=NA
          skyRMS=NA
          skyRMSpval=NA
          tempybak<-NA
          tempxbak<-NA
          templtybak<-NA
          tempylimsbak<-NA
          skybak<-NA
          if (run==1) { medianStat<-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
          if (run==2) { meanStat  <-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
        }
      }
      lo<-min(tempval,na.rm=TRUE)
      hi<-max(tempval,na.rm=TRUE)
      if (!(is.finite(lo)&is.finite(hi))) {
        lo=-5
        hi=5
      }
      xlocs<-floor(pixlocx+(-cuthi:cuthi))
      ylocs<-floor(pixlocy+(-cuthi:cuthi))
      #Select only pixels which are inside the image bounds
      xsel<-which(xlocs>0 & xlocs<=length(origim[,1]))
      ysel<-which(ylocs>0 & ylocs<=length(origim[1,]))
      #Trim to above
      xlocs<-xlocs[xsel]
      ylocs<-ylocs[ysel]
      #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
      if(rem.mask){
        maskim[which(maskim==0)]<-NA
        tempim<-origim*maskim
      }
      #All ref pixels for new image
      tempref<-as.matrix(expand.grid(xlocs,ylocs))
      #Corresponding radii for new pixels from the object of interest
      temprad<-sqrt((tempref[,1]-pixlocx)^2+(tempref[,2]-pixlocy)^2)
      #Keep only pixels inside the radius bounds given by cutlo and cuthi
      keep<-temprad>cutlo & temprad<cuthi
      #Trim
      tempref<-tempref[keep,]
      tempval<-tempim[tempref]
      temprad<-temprad[keep]
      if (length(which(tempval!=0))==0) {
        next
      }
      #Plot
      if (toFile) {
        PlotPNG(file=file.path(path,paste(id,"_skyback.png",sep="")),height=5*res,width=10*res,res=res,pointsize=14)
      }
      layout(cbind(1,2))
      par(mar=c(3.1,3.1,1.1,1.1))
      image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,matrix(magmap(tempim,stretch='asinh',lo=lo,hi=hi,type='num')$map,nrow=nrow(tempim),ncol=ncol(tempim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="",asp=1,useRaster=TRUE,xlim=c(-cuthi,cuthi),ylim=c(-cuthi,cuthi))
      #image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,log10(1-maskim),col=hsv(0,0,0,alpha=1),add=TRUE,useRaster=TRUE)
      points(x=(x.pix-pixloc[1]+1),y=(y.pix-pixloc[2]+1), pch=3)
      if (length(tempxbak)!=0 & !any(is.na(templtybak))) {
        for (bin in 1:length(tempxbak)) { lines(ellipse(a=tempxbak[bin],b=tempxbak[bin],xcen=0,ycen=0),col=grey(0.5),lty=templtybak[bin],lwd=bin.lwd) }
      }
      if (length(tempx)!=0 & !any(is.na(templty))) {
        for (bin in 1:length(tempx)) { lines(ellipse(a=tempx[bin],b=tempx[bin],xcen=0,ycen=0),col='purple',lty=templty[bin],lwd=bin.lwd) }
      }
      magaxis(side=3:4,labels=FALSE)
      magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
      inten<-magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map
      inten[which(is.na(inten))]<-1
      if (length(tempval)>1E4) {
        #Thin out
        samp<-sample(length(tempval),1E4)
        magplot(x=temprad[samp],y=tempval[samp],pch=20,xlab='Radius (pix); Thinned to 1E4 points',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad[samp],na.rm=TRUE)),100,max(temprad[samp],na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.0,1.0)*ifelse(!is.finite(skyRMS),mad(origim),skyRMS),col=hsv(inten[samp]))
        magaxis(side=3:4,labels=FALSE)
      } else {
        magplot(x=temprad,y=tempval,pch=20,xlab='Radius (pix)',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad,na.rm=T)),100,max(temprad,na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.0,1.0)*ifelse(!is.finite(skyRMS),mad(origim),skyRMS),col=hsv(inten))
        magaxis(side=3:4,labels=FALSE)
      }
      if (!any(is.na(templtybak))) {
        abline(v=tempxbak,col=grey(0.5),lty=templtybak,lwd=bin.lwd)
        abline(v=tempx,col='purple',lty=templty,lwd=bin.lwd)
        lines(x=tempxbak,y=tempybak,lwd=all.est.lwd,col=grey(0.5))
        lines(x=tempx,y=tempy,lwd=est.lwd)
        lines(x=tempxbak,y=tempylimsbak[,1],lwd=all.est.lwd,lty=2,col=grey(0.5))
        lines(x=tempxbak,y=tempylimsbak[,2],lwd=all.est.lwd,lty=2,col=grey(0.5))
        lines(x=tempx,y=tempylims[,1],lwd=est.lwd,lty=2)
        lines(x=tempx,y=tempylims[,2],lwd=est.lwd,lty=2)
        abline(h=0,col='darkgreen')
        abline(h=meanStat$sky,col='red')
        abline(h=medianStat$sky,col='red',lty=2)
      }
      dev.off()
    }
  } else {
    for (r in rand) {
      cutlo=cutlovec[r]
      cuthi=cuthivec[r]
      sxl=data.stamp.lims[r,1]
      syl=data.stamp.lims[r,3]
      pixlocx=x.pix[r]
      pixlocy=y.pix[r]
      id=cat.id[r]

      pixloc=c(pixlocx,pixlocy)
      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
      for (run in 2:1) {
        #Find extreme pixels to cut out
        xlocs<-floor(pixlocx+(-cuthi:cuthi))
        ylocs<-floor(pixlocy+(-cuthi:cuthi))
        #Select only pixels which are inside the image bounds
        xsel<-which(xlocs>0 & xlocs<=length(data.stamp[,1]))
        ysel<-which(ylocs>0 & ylocs<=length(data.stamp[1,]))
        #Trim to above
        xlocs<-xlocs[xsel]
        ylocs<-ylocs[ysel]
        #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
        if(rem.mask){
          tempim<-data.stamp*mask.stamp
          tempim[which(mask.stamp==0)]<-NA
        }else{
          tempim<-data.stamp
        }
        #All ref pixels for new image
        tempref<-as.matrix(expand.grid(xlocs,ylocs))
        #Corresponding radii for new pixels from the object of interest
        temprad<-sqrt((tempref[,1]-pixlocx)^2+(tempref[,2]-pixlocy)^2)
        #Keep only pixels inside the radius bounds given by cutlo and cuthi
        keep<-temprad>cutlo & temprad<cuthi
        #Trim
        tempref<-tempref[keep,]
        tempval<-tempim[tempref]
        temprad<-temprad[keep]
        #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
        #if(rem.mask){
        #  temprad<-temprad[tempval!=0]
        #  tempval<-tempval[tempval!=0]
        #}
        #Do iterative <probcut> pixel clipping
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
        templty<-rep(1:2,ceiling(length(tempx)/2))[1:length(tempx)]
        num.in.bin<-tempmedian$Nbins
        #Remove bins with 1 or less skypixels present
        if (any(num.in.bin<=1)) {
          ind<-which(!num.in.bin<=1)
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
            tempybak<-tempy
            tempxbak<-tempx
            templtybak<-templty
            tempylimsbak<-tempylims
            skybak<-sky
            while(any(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr))) & all(!(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))==FALSE){
              nloop<-nloop+1
              tempy<-tempy[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              tempx<-tempx[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              templty<-templty[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              weights<-weights[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)]
              tempylims<-rbind(tempylims[tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr),])
              sky<-sum(tempy*weights)/(sum(weights))
              if (nloop>15) {
                warning("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
                message("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
                break
              }
            }
            #Find the number of running medians that agree with the final sky within error bounds (max=10)
            Nnearsky<-length(which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))
          } else {
            Nnearsky<-1
          }
          if (run==1) { medianStat<-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
          if (run==2) { meanStat  <-data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval) }
        } else {
          #No Sky estimate available - return NAs
          sky=NA
          skyerr=NA
          Nnearsky=NA
          skyRMS=NA
          skyRMSpval=NA
          tempybak<-NA
          tempxbak<-NA
          templtybak<-NA
          tempylimsbak<-NA
          skybak<-sky
          if (run==1) { medianStat<-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
          if (run==2) { meanStat  <-data.frame(sky=NA,skyerr=NA,Nnearsky=NA,skyRMS=NA,skyRMSpval=NA) }
        }
      }
      lo<-min(tempval,na.rm=TRUE)
      hi<-max(tempval,na.rm=TRUE)
      if (!(is.finite(lo)&is.finite(hi))) {
        lo=-5
        hi=5
      }
      xlocs<-floor(pixlocx+(-cuthi:cuthi))
      ylocs<-floor(pixlocy+(-cuthi:cuthi))
      #Select only pixels which are inside the image bounds
      xsel<-which(xlocs>0 & xlocs<=length(data.stamp[,1]))
      ysel<-which(ylocs>0 & ylocs<=length(data.stamp[1,]))
      #Trim to above
      xlocs<-xlocs[xsel]
      ylocs<-ylocs[ysel]
      #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
      if(rem.mask){
        tempim<-data.stamp*mask.stamp
        tempim[which(mask.stamp==0)]<-NA
      }
      #All ref pixels for new image
      tempref<-as.matrix(expand.grid(xlocs,ylocs))
      #Corresponding radii for new pixels from the object of interest
      temprad<-sqrt((tempref[,1]-pixlocx)^2+(tempref[,2]-pixlocy)^2)
      #Keep only pixels inside the radius bounds given by cutlo and cuthi
      keep<-temprad>cutlo & temprad<cuthi
      #Trim
      tempref<-tempref[keep,]
      tempval<-tempim[tempref]
      temprad<-temprad[keep]
      if (length(which(tempval!=0))==0) {
        next
      }
      #Plot
      if (toFile) {
        PlotPNG(file=file.path(path,paste(id,"_skyback.png",sep="")),height=5*res,width=10*res,res=res,pointsize=14)
      }
      layout(cbind(1,2))
      par(mar=c(3.1,3.1,1.1,1.1))
      image(x=1:length(data.stamp[,1])-pixlocx,y=1:length(data.stamp[1,])-pixlocy,matrix(magmap(tempim,stretch='asinh',lo=lo,hi=hi,type='num')$map,nrow=nrow(tempim),ncol=ncol(tempim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="",asp=1,useRaster=TRUE, xlim=c(-cuthi,cuthi),ylim=c(-cuthi,cuthi))
      points(x=(x.pix-pixloc[1]+1),y=(y.pix-pixloc[2]+1), pch=3)
      #if (remmask) {
      #  image(x=1:length(data.stamp[,1])-pixlocx,y=1:length(data.stamp[1,])-pixlocy,log10(1-mask.stamp),col=hsv(0,0,0,alpha=1),add=TRUE,useRaster=TRUE)
      #}
      if (length(tempxbak)!=0 & !any(is.na(templtybak))) {
        for (bin in 1:length(tempxbak)) { lines(ellipse(a=tempxbak[bin],b=tempxbak[bin],xcen=0,ycen=0),col=grey(0.5),lty=templtybak[bin],lwd=bin.lwd) }
      }
      if (length(tempx)!=0 & !any(is.na(templty))) {
        for (bin in 1:length(tempx)) { lines(ellipse(a=tempx[bin],b=tempx[bin],xcen=0,ycen=0),col='purple',lty=templty[bin],lwd=bin.lwd) }
      }
      magaxis(side=3:4,labels=FALSE)
      magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
      inten<-magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map
      inten[which(is.na(inten))]<-1
      if (length(tempval)>1E4) {
        #Thin out
        samp<-sample(length(tempval),1E4)
        magplot(x=temprad[samp],y=tempval[samp],pch=20,xlab='Radius (pix); Thinned to 1E4 points',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad[samp],na.rm=TRUE)),100,max(temprad[samp],na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.0,1.0)*ifelse(!is.finite(skyRMS),mad(data.stamp),skyRMS),col=hsv(inten[samp]))
        magaxis(side=3:4,labels=FALSE)
      } else {
        magplot(x=temprad,y=tempval,pch=20,xlab='Radius (pix)',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad,na.rm=TRUE)),100,max(temprad,na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.0,1.0)*ifelse(!is.finite(skyRMS),mad(data.stamp),skyRMS),col=hsv(inten))
        magaxis(side=3:4,labels=FALSE)
      }
      if (!any(is.na(templtybak))) {
        abline(v=tempxbak,col=grey(0.5),lty=templtybak,lwd=bin.lwd)
        abline(v=tempx,col='purple',lty=templty,lwd=bin.lwd)
        lines(x=tempxbak,y=tempybak,lwd=all.est.lwd,col=grey(0.5))
        lines(x=tempx,y=tempy,lwd=est.lwd)
        lines(x=tempxbak,y=tempylimsbak[,1],lwd=all.est.lwd,lty=2,col=grey(0.5))
        lines(x=tempxbak,y=tempylimsbak[,2],lwd=all.est.lwd,lty=2,col=grey(0.5))
        lines(x=tempx,y=tempylims[,1],lwd=est.lwd,lty=2)
        lines(x=tempx,y=tempylims[,2],lwd=est.lwd,lty=2)
        abline(h=0,col='darkgreen')
        abline(h=meanStat$sky,col='red')
        abline(h=medianStat$sky,col='red',lty=2)
      }
      if (toFile) { dev.off() }
    }
  }
  #Output the foreach data
  return=data.frame(sky=medianStat$sky,skyerr=medianStat$skyerr,Nnearsky=medianStat$Nnearsky,skyRMS=medianStat$skyRMS,skyRMSpval=medianStat$skyRMSpval,
                    sky.mean=meanStat$sky,skyerr.mean=meanStat$skyerr,Nnearsky.mean=meanStat$Nnearsky,skyRMS.mean=meanStat$skyRMS,skyRMSpval.mean=meanStat$skyRMSpval)
}

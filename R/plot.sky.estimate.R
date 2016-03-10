plot.sky.estimate<-function(cat.id,x.pix,y.pix,stamplims=NULL,cutlo=0,cuthi=100,data.stamp,mask.stamp=NULL,remmask=TRUE,radweight=1,clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,probcut=3,plot.all=FALSE,path=NULL,toFile=TRUE,res=120){
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
  if(is.null(stamplims)){
    if (is.matrix(data.stamp)) {
      stamplims<-matrix(c(1,length(data.stamp[,1]),1,length(data.stamp[1,])),ncol=4,nrow=length(x.pix))
    } else {
      stop('Stamplims not provided!')
    }
  }
  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(probcut)
  if (toFile) {
    dir.create(file.path(path,"SkyBackIms"),showWarnings=FALSE)
    path=file.path(path,"SkyBackIms")
  }

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
      sxl=stamplims[r,1]
      syl=stamplims[r,3]
      pixlocx=x.pix[r]
      pixlocy=y.pix[r]
      id=cat.id[r]

      pixloc=c(pixlocx,pixlocy)
      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
      for (run in 1:2) {
        #Find extreme pixels to cut out
        xlocs<-floor(pixlocx+(-cuthi:cuthi))
        ylocs<-floor(pixlocy+(-cuthi:cuthi))
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
        templty<-rep(1:2,ceiling(length(tempx)/2))[1:length(tempx)]
        #Remove bins with 1 or less skypixels present
        if (any(is.na(tempylims)|(tempylims[,2]==tempylims[,1]))) {
          ind<-which((!is.na(tempylims[,1]))&(tempylims[,2]!=tempylims[,1]))
          tempy  <-tempy[ind]
          tempx  <-tempx[ind]
          templty<-templty[ind]
          tempylims<-matrix(tempylims[which((!is.na(tempylims))&(tempylims[,2]!=tempylims[,1]),arr.ind=TRUE)], ncol=2)
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
            if (nloop==1E8) {
              warning("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
              message("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
              break
            }
          }
          #Find the number of running medians that agree with the final sky within error bounds (max=10)
          Nnearsky<-length(which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))
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
      if(remmask){
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
        CairoPNG(file=file.path(path,paste(id,"_skyback.png",sep="")),height=6*res,width=10*res,res=res)
      }
      layout(cbind(1,2))
      image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,magmap(origim,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="",asp=1,useRaster=TRUE,xlim=c(-cuthi,cuthi),ylim=c(-cuthi,cuthi))
      image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,log10(1-maskim),col=hsv(0,0,0,alpha=1),add=TRUE,useRaster=TRUE)
      if (length(tempxbak)!=0 & !any(is.na(templtybak))) {
        for (bin in 1:length(tempxbak)) { lines(ellipse(a=tempxbak[bin],b=tempxbak[bin],xcen=0,ycen=0),col='darkgrey',lty=templtybak[bin],lwd=3) }
      }
      if (length(tempx)!=0 & !any(is.na(templty))) {
        for (bin in 1:length(tempx)) { lines(ellipse(a=tempx[bin],b=tempx[bin],xcen=0,ycen=0),col='purple',lty=templty[bin],lwd=3) }
      }
      magaxis(side=1:4,labels=FALSE)
      magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
      points(x=(x.pix-pixloc[1]+1),y=(y.pix-pixloc[2]+1), pch=3)
      inten<-magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map
      inten[which(is.na(inten))]<-1
      if (length(tempval)>1E4) {
        #Thin out
        samp<-sample(length(tempval),1E4)
        magplot(x=temprad[samp],y=tempval[samp],pch=20,xlab='Radius (pix); Thinned to 1E4 points',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad[samp],na.rm=TRUE)),100,max(temprad[samp],na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.5,1.5)*ifelse(!is.finite(skyRMS),sd(origim),skyRMS),col=hsv(inten[samp]))
      } else {
        magplot(x=temprad,y=tempval,pch=20,xlab='Radius (pix)',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad,na.rm=T)),100,max(temprad,na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.5,1.5)*ifelse(!is.finite(skyRMS),sd(origim),skyRMS),col=hsv(inten))
      }
      abline(v=tempxbak,col='darkgrey',lty=templtybak,lwd=3)
      abline(v=tempx,col='purple',lty=templty,lwd=3)
      lines(x=tempxbak,y=tempybak,lwd=1,col='darkgrey')
      lines(x=tempx,y=tempy,lwd=2)
      lines(x=tempxbak,y=tempylimsbak[,1],lwd=1,lty=2,col='darkgrey')
      lines(x=tempxbak,y=tempylimsbak[,2],lwd=1,lty=2,col='darkgrey')
      lines(x=tempx,y=tempylims[,1],lwd=2,lty=2)
      lines(x=tempx,y=tempylims[,2],lwd=2,lty=2)
      abline(h=0,col='darkgreen')
      abline(h=meanStat$sky,col='red')
      abline(h=medianStat$sky,col='red',lty=4)
      dev.off()
    }
  } else {
    for (r in rand) {
      cutlo=cutlovec[r]
      cuthi=cuthivec[r]
      sxl=stamplims[r,1]
      syl=stamplims[r,3]
      pixlocx=x.pix[r]
      pixlocy=y.pix[r]
      id=cat.id[r]

      pixloc=c(pixlocx,pixlocy)
      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
      for (run in 1:2) {
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
        if(remmask){
          tempim<-data.stamp*mask.stamp
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
        templty<-rep(1:2,ceiling(length(tempx)/2))[1:length(tempx)]
        #Remove bins with 1 or less skypixels present
        if (any(is.na(tempylims)|(tempylims[,2]==tempylims[,1]))) {
          ind<-which((!is.na(tempylims[,1]))&(tempylims[,2]!=tempylims[,1]))
          tempy  <-tempy[ind]
          tempx  <-tempx[ind]
          templty<-templty[ind]
          tempylims<-matrix(tempylims[which((!is.na(tempylims))&(tempylims[,2]!=tempylims[,1]),arr.ind=TRUE)], ncol=2)
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
            if (nloop==1E8) {
              warning("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
              message("No running medians are inside the 1-sigma bound of the sky after 1E8 loops. Breaking.")
              break
            }
          }
          #Find the number of running medians that agree with the final sky within error bounds (max=10)
          Nnearsky<-length(which(tempylims[,1]<=(sky+skyerr) & tempylims[,2]>=(sky-skyerr)))
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
      if(remmask){
        tempim<-data.stamp*mask.stamp
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
        CairoPNG(file=file.path(path,paste(id,"_skyback.png",sep="")),height=6*res,width=10*res,res=res)
      }
      layout(cbind(1,2))
      image(x=1:length(data.stamp[,1])-pixlocx,y=1:length(data.stamp[1,])-pixlocy,magmap(data.stamp,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="",asp=1,useRaster=TRUE, xlim=c(-cuthi,cuthi),ylim=c(-cuthi,cuthi))
      if (remmask) {
        image(x=1:length(data.stamp[,1])-pixlocx,y=1:length(data.stamp[1,])-pixlocy,log10(1-mask.stamp),col=hsv(0,0,0,alpha=1),add=TRUE,useRaster=TRUE)
      }
      if (length(tempxbak)!=0 & !any(is.na(templtybak))) {
        for (bin in 1:length(tempxbak)) { lines(ellipse(a=tempxbak[bin],b=tempxbak[bin],xcen=0,ycen=0),col='darkgrey',lty=templtybak[bin],lwd=3) }
      }
      if (length(tempx)!=0 & !any(is.na(templty))) {
        for (bin in 1:length(tempx)) { lines(ellipse(a=tempx[bin],b=tempx[bin],xcen=0,ycen=0),col='purple',lty=templty[bin],lwd=3) }
      }
      magaxis(side=1:4,labels=FALSE)
      magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
      points(x=(x.pix-pixloc[1]+1),y=(y.pix-pixloc[2]+1), pch=3)
      inten<-magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map
      inten[which(is.na(inten))]<-1
      if (length(tempval)>1E4) {
        #Thin out
        samp<-sample(length(tempval),1E4)
        magplot(x=temprad[samp],y=tempval[samp],pch=20,xlab='Radius (pix); Thinned to 1E4 points',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad[samp],na.rm=TRUE)),100,max(temprad[samp],na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.5,1.5)*ifelse(!is.finite(skyRMS),sd(origim),skyRMS),col=hsv(inten[samp]))
      } else {
        magplot(x=temprad,y=tempval,pch=20,xlab='Radius (pix)',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad,na.rm=TRUE)),100,max(temprad,na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.5,1.5)*ifelse(!is.finite(skyRMS),sd(origim),skyRMS),col=hsv(inten))
      }
      if (!any(is.na(templtybak))) {
        abline(v=tempxbak,col='darkgrey',lty=templtybak,lwd=3)
        abline(v=tempx,col='purple',lty=templty,lwd=3)
        lines(x=tempxbak,y=tempybak,lwd=1,col='darkgrey')
        lines(x=tempx,y=tempy,lwd=2)
        lines(x=tempxbak,y=tempylimsbak[,1],lwd=1,lty=2,col='darkgrey')
        lines(x=tempxbak,y=tempylimsbak[,2],lwd=1,lty=2,col='darkgrey')
        lines(x=tempx,y=tempylims[,1],lwd=2,lty=2)
        lines(x=tempx,y=tempylims[,2],lwd=2,lty=2)
        abline(h=0,col='darkgreen')
        abline(h=meanStat$sky,col='red')
        abline(h=medianStat$sky,col='red',lty=4)
      }
      if (toFile) { dev.off() }
    }
  }
  #Output the foreach data
  return=NULL
}

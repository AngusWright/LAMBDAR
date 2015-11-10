PlotSkyback<-function(id_g,x_p,y_p,stamplims,cutlo=0,cuthi=100,im_mask,imm_mask,remmask=TRUE,radweight=1,clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,probcut=3,plotall=FALSE,path=NULL,mpiopts=mpiopts){
  cutup=TRUE
  if(length(x_p) != length(y_p)){stop('x_pix and y_pix lengths do not match!')}
  if(length(x_p) != length(cutlo)){stop('x_pix and cutlo lengths do not match!')}
  if(length(x_p) != length(cuthi)){stop('x_pix and cuthi lengths do not match!')}
  if(length(x_p) != length(im_mask)){
    if (is.matrix(im_mask)) {
      cutup=FALSE
    } else {
      stop('x_pix and im_mask lengths do not match!')
    }
  }
  cutlo[cutlo<hardlo*PSFFWHMinPIX]<-hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]<-hardhi*PSFFWHMinPIX
  cutlovec<-round(cutlo)
  cuthivec<-round(cuthi)
  probcut<-1-pnorm(probcut)
  dir.create(file.path(path,"SkyBackIms"),showWarnings=FALSE)
  path=file.path(path,"SkyBackIms")

  if(plotall) {
    rand<-1:length(x_p)
  } else {
    rand<-sample(length(x_p),min(10,length(x_p)))
  }
  if (cutup) {
    #output<-foreach(id=id_g, cutlo=cutlovec, cuthi=cuthivec, pixlocx=x_p, pixlocy=y_p,sxl=stamplims[,1],syl=stamplims[,3], origim=im_mask, maskim=imm_mask, .export=c('x_p','y_p'), .options.mpi=mpiopts, .combine='rbind') %dopar% {
    for (r in rand) {
      cutlo=cutlovec[r]
      cuthi=cuthivec[r]
      origim=im_mask[[r]]
      maskim=imm_mask[[r]]
      sxl=stamplims[r,1]
      syl=stamplims[r,3]
      pixlocx=x_p[r]
      pixlocy=y_p[r]
      id=id_g[r]

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
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T,type='mean') }
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
      lo<-min(tempval,na.rm=T)
      hi<-max(tempval,na.rm=T)
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
      #useCairo=TRUE
      #if (useCairo) {
      #  #if (r==rand[1]) { cat("..using Cairo PDF...") }
      #  library("Cairo",lib.loc='/gamah/awright/src/R/lib/')
      #  CairoPDF(file=file.path(path,paste(id,"_skyback.pdf",sep="")),height=7,width=10)
      #} else {
        pdf(file=file.path(path,paste(id,"_skyback.pdf",sep="")),height=7,width=10)
      #}
      layout(cbind(1,2))
      #image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,magmap(origim,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(0,0,seq(2/3,0,length=256)),axes=F,ylab="",xlab="",main="",asp=1,useRaster=TRUE)
      image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,magmap(origim,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(seq(2/3,0,length=256)),axes=F,ylab="",xlab="",main="",asp=1,useRaster=TRUE,xlim=c(-cuthi,cuthi),ylim=c(-cuthi,cuthi))
      image(x=1:length(origim[,1])-pixlocx,y=1:length(origim[1,])-pixlocy,log10(1-maskim),col=hsv(0,0,0,alpha=1),add=TRUE,useRaster=TRUE)
      if (length(tempxbak)!=0 & !any(is.na(templtybak))) {
        for (bin in 1:length(tempxbak)) { lines(ellipse(a=tempxbak[bin],b=tempxbak[bin],xcen=0,ycen=0),col='darkgrey',lty=templtybak[bin],lwd=3) }
      }
      if (length(tempx)!=0 & !any(is.na(templty))) {
        for (bin in 1:length(tempx)) { lines(ellipse(a=tempx[bin],b=tempx[bin],xcen=0,ycen=0),col='purple',lty=templty[bin],lwd=3) }
      }
      magaxis(side=1:4,labels=F)
      magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
      points(x=(x_p-pixloc[1]+1),y=(y_p-pixloc[2]+1), pch=3)
      inten<-magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map
      inten[which(is.na(inten))]<-1
      magplot(x=temprad,y=tempval,pch=20,xlab='Radius (pix)',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad,na.rm=T)),100,max(temprad,na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.5,1.5)*ifelse(!is.finite(skyRMS),sd(origim),skyRMS),col=hsv(inten))
      #points(x=temprad,y=tempval,pch=20,col=hsv(magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map))
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
    #output<-foreach(id=id_g, cutlo=cutlovec, cuthi=cuthivec, pixlocx=x_p, pixlocy=y_p, sxl=stamplims[,1],syl=stamplims[,3],.export=c('x_p','y_p'), .options.mpi=mpiopts, .combine='rbind') %dopar% {
    for (r in rand) {
      cutlo=cutlovec[r]
      cuthi=cuthivec[r]
      sxl=stamplims[r,1]
      syl=stamplims[r,3]
      pixlocx=x_p[r]
      pixlocy=y_p[r]
      id=id_g[r]

      pixloc=c(pixlocx,pixlocy)
      pixlocx=pixlocx-(sxl-1)
      pixlocy=pixlocy-(syl-1)
      for (run in 1:2) {
        #Find extreme pixels to cut out
        xlocs<-floor(pixlocx+(-cuthi:cuthi))
        ylocs<-floor(pixlocy+(-cuthi:cuthi))
        #Select only pixels which are inside the image bounds
        xsel<-which(xlocs>0 & xlocs<=length(im_mask[,1]))
        ysel<-which(ylocs>0 & ylocs<=length(im_mask[1,]))
        #Trim to above
        xlocs<-xlocs[xsel]
        ylocs<-ylocs[ysel]
        #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
        if(remmask){
          tempim<-im_mask*imm_mask
        }else{
          tempim<-im_mask
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
        if (run==1) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T) }
        if (run==2) { tempmedian<-magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T,type='mean') }
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
      lo<-min(tempval,na.rm=T)
      hi<-max(tempval,na.rm=T)
      if (!(is.finite(lo)&is.finite(hi))) {
        lo=-5
        hi=5
      }
      xlocs<-floor(pixlocx+(-cuthi:cuthi))
      ylocs<-floor(pixlocy+(-cuthi:cuthi))
      #Select only pixels which are inside the image bounds
      xsel<-which(xlocs>0 & xlocs<=length(im_mask[,1]))
      ysel<-which(ylocs>0 & ylocs<=length(im_mask[1,]))
      #Trim to above
      xlocs<-xlocs[xsel]
      ylocs<-ylocs[ysel]
      #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
      if(remmask){
        tempim<-im_mask*imm_mask
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
      #useCairo=TRUE
      #if (useCairo) {
      #  #if (r==rand[1]) { cat("..using Cairo PDF...") }
      #  library("Cairo",lib.loc='/gamah/awright/src/R/lib/')
      #  CairoPDF(file=file.path(path,paste(id,"_skyback.pdf",sep="")),height=7,width=10)
      #} else {
        pdf(file=file.path(path,paste(id,"_skyback.pdf",sep="")),height=7,width=10)
      #}
      layout(cbind(1,2))
      #image(x=1:length(im_mask[,1])-pixlocx,y=1:length(im_mask[1,])-pixlocy,magmap(im_mask,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(0,0,seq(0,1,length=256)),axes=F,ylab="",xlab="",main="",asp=1,useRaster=TRUE)
      #image(x=1:length(im_mask[,1])-pixlocx,y=1:length(im_mask[1,])-pixlocy,magmap(tempim,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(seq(2/3,0,length=256)),axes=F,ylab="",xlab="",main="",asp=1,useRaster=TRUE,add=TRUE)
      image(x=1:length(im_mask[,1])-pixlocx,y=1:length(im_mask[1,])-pixlocy,magmap(im_mask,stretch='asinh',lo=lo,hi=hi,type='num')$map,col=hsv(seq(2/3,0,length=256)),axes=F,ylab="",xlab="",main="",asp=1,useRaster=TRUE, xlim=c(-cuthi,cuthi),ylim=c(-cuthi,cuthi))
      image(x=1:length(im_mask[,1])-pixlocx,y=1:length(im_mask[1,])-pixlocy,log10(1-imm_mask),col=hsv(0,0,0,alpha=1),add=TRUE,useRaster=TRUE)
      if (length(tempxbak)!=0 & !any(is.na(templtybak))) {
        for (bin in 1:length(tempxbak)) { lines(ellipse(a=tempxbak[bin],b=tempxbak[bin],xcen=0,ycen=0),col='darkgrey',lty=templtybak[bin],lwd=3) }
      }
      if (length(tempx)!=0 & !any(is.na(templty))) {
        for (bin in 1:length(tempx)) { lines(ellipse(a=tempx[bin],b=tempx[bin],xcen=0,ycen=0),col='purple',lty=templty[bin],lwd=3) }
      }
      magaxis(side=1:4,labels=F)
      magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
      points(x=(x_p-pixloc[1]+1),y=(y_p-pixloc[2]+1), pch=3)
      inten<-magmap(tempval,lo=lo,hi=hi,type='num',range=c(0,2/3),flip=TRUE,stretch='asinh')$map
      inten[which(is.na(inten))]<-1
      magplot(x=temprad,y=tempval,pch=20,xlab='Radius (pix)',ylab="Pixel Value",xlim=c(0,ifelse(!is.finite(max(temprad,na.rm=T)),100,max(temprad,na.rm=TRUE))),ylim=ifelse(!is.finite(sky),0,sky)+c(-1.5,1.5)*ifelse(!is.finite(skyRMS),sd(im_mask),skyRMS),col=hsv(inten))
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
      dev.off()
    }
  }
  #Output the foreach data
  return=NULL
}

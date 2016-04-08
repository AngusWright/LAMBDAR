plot.ran.cor<-function(data.stamp,ap.stamp,mask.stamp=NULL,ap.stamp.lims=NULL,data.stamp.lims=NULL,mask.stamp.lims=NULL,toFile=FALSE,rem.mask=FALSE,numIters=1E2,path="./",plot.all=FALSE,sigclip=3,nclip=3,res=120,cat.id=NULL,cat.x=NULL,cat.y=NULL){
  if(is.matrix(ap.stamp)) {
    #We have 1 object
    ap.stamp<-list(ap.stamp)
  }
  #Check Ap Stamp limits
  if (is.null(ap.stamp.lims)){
    #Ap stamps must be the same dimension as the image!
    warning('ap.stamp.lims is not provided, so we assume the limits are the stamp edges')
    ap.stamp.lims<-rbind(foreach(ap=ap.stamp,.combine='rbind')%dopar%{return=c(1,length(ap[,1]),1,length(ap[1,]))})
  } else if (length(ap.stamp) != length(ap.stamp.lims[,1])){
    stop('ap.stamp.lims is not the same length as ap.stamp!')
  }
  #Check Data stamp Limits
  if (is.null(data.stamp.lims)){
    #Data stamps must be the same dimension as the image!
    warning('data.stamp.lims is not provided, so we assume the limits are the stamp edges')
    if (is.matrix(data.stamp)) {
      data.stamp.lims<-matrix(c(1,length(data.stamp[,1]),1,length(data.stamp[1,])),nrow=length(ap.stamp),ncol=4,byrow=T)
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
    #Check mask stamp Limits
    if (is.null(mask.stamp.lims)){
      #Mask stamps must be the same dimension as the image!
      warning('mask.stamp.lims is not provided, so we assume the limits are the stamp edges')
      if (is.matrix(mask.stamp)) {
        mask.stamp.lims<-matrix(c(1,length(mask.stamp[,1]),1,length(mask.stamp[1,])),nrow=length(ap.stamp),ncol=4,byrow=T)
      } else {
        mask.stamp.lims<-rbind(foreach(ap=mask.stamp,.combine='rbind')%dopar%{return=c(1,length(ap[,1]),1,length(ap[1,]))})
      }
    } else if (length(mask.stamp) != length(mask.stamp.lims[,1])){
      stop('mask.stamp.lims is not the same length as mask.stamp!')
    }
  }
  if (rem.mask) {
    if(length(ap.stamp) != length(mask.stamp.lims[,1])){
      stop('ap.stamp and mask.stamp.lims lengths do not match!')
    }
  }
  if(length(ap.stamp) != length(data.stamp.lims[,1])){
    stop('ap.stamp and data.stamp.lims lengths do not match!')
  }
  cutup<-TRUE
  if(length(ap.stamp) != length(data.stamp)){
    if (is.matrix(data.stamp)) {
      cutup<-FALSE
    } else {
      stop("ap.stamp and data.stamp lengths do not match!")
    }
  }
  if (rem.mask) {
    if(length(ap.stamp) != length(mask.stamp)) {
      if (cutup){
        stop("ap.stamp and mask.stamp lengths do not match!")
      } else if (!is.matrix(mask.stamp)) {
        stop("data.stamp and mask.stamp are not of the same form (both matrix or both list)")
      }
    }
  }
  if (toFile) {
    if (is.null(cat.id)) { stop("Output to file requires cat.id to be specified (used for filenames)") }
    if (!rem.mask) {
      dir.create(file.path(path,"ran.corIms"),showWarnings=FALSE)
      path=file.path(path,"ran.corIms")
    } else {
      dir.create(file.path(path,"BlanksCorIms"),showWarnings=FALSE)
      path=file.path(path,"BlanksCorIms")
    }
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

  if(plot.all) {
    rand<-1:length(x.pix)
  } else {
    rand<-sample(length(x.pix),min(10,length(x.pix)))
  }
  if (cutup) {
    for (i in rand) {
     origim=data.stamp[[i]]
    maskim=mask.stamp[[i]]
         ap=ap.stamp[[i]]
     sxl=ap.stamp.lims[i,1]
     sxh=ap.stamp.lims[i,2]
     syl=ap.stamp.lims[i,3]
     syh=ap.stamp.lims[i,4]
     imsxl=data.stamp.lims[i,1]
     imsxh=data.stamp.lims[i,2]
     imsyl=data.stamp.lims[i,3]
     imsyh=data.stamp.lims[i,4]

      #Get number of cols & rows in image stamp
      nc<-ncol(origim)
      nr<-nrow(origim)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      #Mask object aperture
      origim[sxl:sxh,syl:syh][which(zapsmall(ap)!=0,arr.ind=TRUE)]<-NA
      #Get Aperture details
      gdap<-which(ap!=0)
      gdapv<-ap[which(ap!=0)]
      #If Remasking
      if(rem.mask){
        #Set mask 0s to NA
        maskim[which(maskim==0)]<-NA
        #Multiply Image and Mask together
        tempim<-origim*maskim
        #For Niters, calculate random Flux
        ranaps<-tempim*0
        tempvec<-matrix(NA,ncol=length(which(ap!=0)),nrow=numIters)
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          ranaps[as.matrix(ind)]<-ranaps[as.matrix(ind)]+gdapv*0.1
          tempvec[iter,]<-tempim[as.matrix(ind)]*gdapv
          #Sum shifted image and Aperture to return Flux
          val<-which(!is.na(tempim[as.matrix(ind)]))
          temp<-tempvec[iter,val]
          val<-gdapv[val]
          sumap[iter]<-sum(val)
          #Sum Flux, correcting for masked pixels
          flux[iter]<-sum(temp*val,na.rm=TRUE)#/(sum(gdapv)/sumap[iter])
        }
        #Return Result
        s<-which(is.finite(flux)&is.finite(sumap))
        if (length(s)!=length(flux)) {
          flux<-flux[s]
          sumap<-sumap[s]
        }
        wflux<-sum(flux*sumap)/sum(sumap)
        wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
        wmad<-weightedMad(flux,sumap)
        if (sigclip>0) {
          for (iter in 1:nclip) {
            ind<-which((flux-wflux)/wsd <= sigclip)
            wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
            wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
            wmad<-weightedMad(flux[ind],sumap[ind])
          }
        }
        dat=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
        lim<-(ceiling(log10(max(abs(quantile(origim,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        if (toFile) { CairoPNG(file=file.path(path,paste(cat.id[i],"_blankscor.png",sep="")),height=6*res,width=10*res,res=res) }
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(origim[,1])-(x.pix[i]-imsxl),y=1:length(origim[1,])-(y.pix[i]-imsyl),matrix(magmap(tempim,lo=0.01,hi=0.99,stretch='asinh')$map,ncol=ncol(tempim),nrow=nrow(tempim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(origim[,1])-(x.pix[i]-imsxl),y=1:length(origim[1,])-(y.pix[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x.pix-x.pix[i]+1),y=(y.pix-y.pix[i]+1), pch=3)
        #Convert pix values onto asinh scale
        #X-Axis Limits; focus on 0 with ± limits around there
        #Stretch data to a useful scale
        stretchscale=ifelse(lim<0,1.5*10^(abs(lim)+1),1.5*10^(-1*(lim-1)))
        #Apply transformation
        tempvecstretch<-magmap(tempvec,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        #Calculate Axes Major Tick Mark Labels
        axespoints<-c(-10^(lim:(lim-2)),0,10^((lim-2):lim))
        #Calculate Axes Major Tick Mark Lenghts
        asinhtcls<-c(rep(0.5,length(axespoints)),rep(0.2,(length(axespoints)-1)*4))
        #Add Minor tick marks
        axespoints<-c(axespoints,c(-10^(lim:(lim-2)),10^((lim-2):lim))/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*2/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*3/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*4/5)
        #Calculate Tick locations on transformed axes
        asinhticks=magmap(axespoints,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',stretchscale=stretchscale)$map
        #Remove overlapping Tick Marks
        if (any(duplicated(asinhticks))) {
          while(length(which(duplicated(asinhticks)))>0) {
            kk<-which(duplicated(asinhticks))[1]
            indall<-which(asinhticks==asinhticks[kk])
            ind2<-indall[which.min(abs(axespoints[indall]))]
            axespoints<-c(axespoints[ind2],axespoints[-1*indall])
            asinhticks<-c(asinhticks[ind2],asinhticks[-1*indall])
            asinhtcls<-c(asinhtcls[ind2],asinhtcls[-1*indall])
          }
        }
        #Get pixel histogram on transformed axes
        pix<-hist(as.numeric(tempvecstretch),plot=FALSE,breaks=seq(0,1,length=100))
        #Plot Histogram and Count axis
        magplot(x=rev(rev(pix$breaks)[-1]),y=pix$counts,type='s',xlab='Randoms Pixel Values (pix)',ylab="Count",side=2,xlim=c(0,1),main="Pixel Histogram",log='y',ylim=c(1,10^(log10(max(pix$counts))+1)))
        #Convert Labels to Pretty style
        labs<-floor(log10(abs(axespoints)))
        pref<-ifelse(axespoints<0,"-","")
        labs<-paste(pref,ifelse(labs>0,"1e+","1e"),labs,"",sep="")
        labs[which(labs=="1e+Inf")]<-"0"
        labs[which(labs=="1e-Inf")]<-"0"
        check = grep("1e+", labs, fixed = TRUE)
        labs[check] = paste(sub("1e+", "10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("1e-", labs, fixed = TRUE)
        labs[check] = paste(sub("1e-", "10^{-", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e+", labs, fixed = TRUE)
        labs[check] = paste(sub("e+", "*x*10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e-", labs, fixed = TRUE)
        labs[check] = paste(sub("e-", "*x*10^{-", labs[check], fixed = TRUE), "}", sep = "")
        axespoints<-parse(text=labs)
        #Draw X major and minor ticks
        axis(1,asinhticks[which(asinhtcls==0.5)],labels=FALSE,tcl=0.5)
        axis(1,asinhticks[which(asinhtcls==0.2)],labels=FALSE,tcl=0.2)
        mtext('Blanks Pixel Values (pix)', 1, line = 2)
        #Draw major tick labels
        ind<-which(asinhtcls==0.5 & (labs=="0" | abs(asinhticks-0.5)>0.1))
        axis(1,asinhticks[ind],labels=axespoints[ind],tcl=0)
        #Draw histograms for each bin.
        for(iter in 1:numIters) {
          if (length(which(!is.na(tempvec[iter,])))>0) {
            tempvecstretch<-magmap(tempvec[iter,],lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
            tmp<-hist(tempvecstretch,plot=FALSE,breaks=pix$breaks)
            lines(x=rev(rev(pix$breaks)[-1]),y=tmp$counts,type='s',col=hsv(seq(2/3,0,length=numIters))[iter])
          }
        }
        abline(v=magmap(dat$randMean.mean/(max(sumap,na.rm=T)),lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Blanks Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",signif(dat$randMean.mean/max(sumap,na.rm=TRUE),digits=3)," (per pix)\nStd Dev = ",signif(dat$randMean.SD/sqrt(max(sumap,na.rm=TRUE)),digits=3)," (per pix)",sep=""),cex=0.6)
        boxx<-magmap(flux/sumap,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.2),boxwex=2)
        if (toFile) { dev.off() }
      }else{
        ranaps<-origim*0
        tempvec<-matrix(NA,ncol=length(which(ap!=0)),nrow=numIters)
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          ranaps[as.matrix(ind)]<-ranaps[as.matrix(ind)]+gdapv*0.1
          tempvec[iter,]<-origim[as.matrix(ind)]*gdapv
          val<-which(!is.na(origim[as.matrix(ind)]))
          temp<-tempvec[iter,val]
          val<-gdapv[val]
          sumap[iter]<-sum(val)
          #Sum Flux, correcting for masked pixels
          flux[iter]<-sum(temp*val,na.rm=TRUE)
        }
        #Return Result
        s<-which(is.finite(flux)&is.finite(sumap))
        if (length(s)!=length(flux)) {
          flux<-flux[s]
          sumap<-sumap[s]
        }
        wflux<-sum(flux*sumap)/sum(sumap)
        wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
        wmad<-weightedMad(flux,sumap)
        if (sigclip>0) {
          for (iter in 1:nclip) {
            ind<-which((flux-wflux)/wsd <= sigclip)
            wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
            wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
            wmad<-weightedMad(flux[ind],sumap[ind])
          }
        }
        dat=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
        lim<-(ceiling(log10(max(abs(quantile(origim,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        if (toFile) { CairoPNG(file=file.path(path,paste(cat.id[i],"_rancor.png",sep="")),height=6*res,width=10*res,res=res) }
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(origim[,1])-(x.pix[i]-imsxl),y=1:length(origim[1,])-(y.pix[i]-imsyl),matrix(magmap(origim,lo=0.01,hi=0.99,stretch='asinh')$map,ncol=ncol(origim),nrow=nrow(origim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(origim[,1])-(x.pix[i]-imsxl),y=1:length(origim[1,])-(y.pix[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x.pix-x.pix[i]+1),y=(y.pix-y.pix[i]+1), pch=3)
        #Convert pix values onto asinh scale
        #X-Axis Limits; focus on 0 with ± limits around there
        #Stretch data to a useful scale
        stretchscale=ifelse(lim<0,1.5*10^(abs(lim)+1),1.5*10^(-1*(lim-1)))
        #Apply transformation
        tempvecstretch<-magmap(tempvec,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        #Calculate Axes Major Tick Mark Labels
        axespoints<-c(-10^(lim:(lim-2)),0,10^((lim-2):lim))
        #Calculate Axes Major Tick Mark Lenghts
        asinhtcls<-c(rep(0.5,length(axespoints)),rep(0.2,(length(axespoints)-1)*4))
        #Add Minor tick marks
        axespoints<-c(axespoints,c(-10^(lim:(lim-2)),10^((lim-2):lim))/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*2/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*3/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*4/5)
        #Calculate Tick locations on transformed axes
        asinhticks=magmap(axespoints,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',stretchscale=stretchscale)$map
        #Remove overlapping Tick Marks
        if (any(duplicated(asinhticks))) {
          while(length(which(duplicated(asinhticks)))>0) {
            kk<-which(duplicated(asinhticks))[1]
            indall<-which(asinhticks==asinhticks[kk])
            ind2<-indall[which.min(abs(axespoints[indall]))]
            axespoints<-c(axespoints[ind2],axespoints[-1*indall])
            asinhticks<-c(asinhticks[ind2],asinhticks[-1*indall])
            asinhtcls<-c(asinhtcls[ind2],asinhtcls[-1*indall])
          }
        }
        #Get pixel histogram on transformed axes
        pix<-hist(as.numeric(tempvecstretch),plot=FALSE,breaks=seq(0,1,length=100))
        #Plot Histogram and Count axis
        magplot(x=rev(rev(pix$breaks)[-1]),y=pix$counts,type='s',xlab='Randoms Pixel Values (pix)',ylab="Count",side=2,xlim=c(0,1),log='y',main="Pixel Histogram",ylim=c(1,10^(log10(max(pix$counts))+1)))
        #Convert Labels to Pretty style
        labs<-floor(log10(abs(axespoints)))
        pref<-ifelse(axespoints<0,"-","")
        labs<-paste(pref,ifelse(labs>0,"1e+","1e"),labs,"",sep="")
        labs[which(labs=="1e+Inf")]<-"0"
        labs[which(labs=="1e-Inf")]<-"0"
        check = grep("1e+", labs, fixed = TRUE)
        labs[check] = paste(sub("1e+", "10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("1e-", labs, fixed = TRUE)
        labs[check] = paste(sub("1e-", "10^{-", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e+", labs, fixed = TRUE)
        labs[check] = paste(sub("e+", "*x*10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e-", labs, fixed = TRUE)
        labs[check] = paste(sub("e-", "*x*10^{-", labs[check], fixed = TRUE), "}", sep = "")
        axespoints<-parse(text=labs)
        #Draw X major and minor ticks
        axis(1,asinhticks[which(asinhtcls==0.5)],labels=FALSE,tcl=0.5)
        axis(1,asinhticks[which(asinhtcls==0.2)],labels=FALSE,tcl=0.2)
        mtext('Randoms Pixel Values (pix)', 1, line = 2)
        #Draw major tick labels
        ind<-which(asinhtcls==0.5 & (labs=="0" | abs(asinhticks-0.5)>0.1))
        axis(1,asinhticks[ind],labels=axespoints[ind],tcl=0)
        #Draw histograms for each bin.
        for(iter in 1:numIters) {
          if (length(which(!is.na(tempvec[iter,])))>0) {
            tempvecstretch<-magmap(tempvec[iter,],lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
            tmp<-hist(tempvecstretch,plot=FALSE,breaks=pix$breaks)
            lines(x=rev(rev(tmp$breaks)[-1]),y=tmp$counts,type='s',col=hsv(seq(2/3,0,length=numIters))[iter])
          }
        }
        abline(v=magmap(dat$randMean.mean/(max(sumap,na.rm=T)),lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Random Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",signif(dat$randMean.mean/max(sumap,na.rm=TRUE),digits=3)," (per pix)\nStd Dev = ",signif(dat$randMean.SD/sqrt(max(sumap,na.rm=TRUE)),digits=3)," (per pix)",sep=""),cex=0.6)
        boxx<-magmap(flux/sumap,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.5),boxwex=2)
        if (toFile) { dev.off() }
      }
    }
  } else {
    for (i in rand) {
         ap=ap.stamp[[i]]
     sxl=ap.stamp.lims[i,1]
     sxh=ap.stamp.lims[i,2]
     syl=ap.stamp.lims[i,3]
     syh=ap.stamp.lims[i,4]
     imsxl=data.stamp.lims[i,1]
     imsxh=data.stamp.lims[i,2]
     imsyl=data.stamp.lims[i,3]
     imsyh=data.stamp.lims[i,4]

      #Get number of cols & rows in image stamp
      nc<-ncol(data.stamp)
      nr<-nrow(data.stamp)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      #Mask object aperture
      tempim<-data.stamp
      tempim[sxl:sxh,syl:syh][which(zapsmall(ap)!=0,arr.ind=TRUE)]<-NA
      #Get Aperture details
      gdap<-which(ap!=0)
      gdapv<-ap[which(ap!=0)]
      #If Remasking
      if(rem.mask){
        #Set mask 0s to NA
        mask.stamp[which(mask.stamp==0)]<-NA
        #Multiply Image and Mask together
        tempim<-tempim*mask.stamp
        #Ftempimers, calculate random Flux
        ranaps<-tempim*0
        tempvec<-matrix(NA,ncol=length(which(ap!=0)),nrow=numIters)
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          ranaps[as.matrix(ind)]<-ranaps[as.matrix(ind)]+gdapv*0.1
          tempvec[iter,]<-tempim[as.matrix(ind)]*gdapv
          val<-which(!is.na(tempim[as.matrix(ind)]))
          temp<-tempvec[iter,val]
          val<-gdapv[val]
          sumap[iter]<-sum(val)
          #Sum Flux, correcting for masked pixels
          flux[iter]<-sum(temp*val,na.rm=TRUE)
        }
        #Return Result
        s<-which(is.finite(flux)&is.finite(sumap))
        if (length(s)!=length(flux)) {
          flux<-flux[s]
          sumap<-sumap[s]
        }
        wflux<-sum(flux*sumap)/sum(sumap)
        wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
        wmad<-weightedMad(flux,sumap)
        if (sigclip>0) {
          for (iter in 1:nclip) {
            ind<-which((flux-wflux)/wsd <= sigclip)
            wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
            wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
            wmad<-weightedMad(flux[ind],sumap[ind])
          }
        }
        dat=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
        lim<-(ceiling(log10(max(abs(quantile(data.stamp,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        if (toFile) { CairoPNG(file=file.path(path,paste(cat.id[i],"_blankscor.png",sep="")),height=6*res,width=10*res,res=res) }
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(data.stamp[,1])-(x.pix[i]-imsxl),y=1:length(data.stamp[1,])-(y.pix[i]-imsyl),matrix(magmap(tempim,lo=0.01,hi=0.99,stretch='asinh')$map,ncol=ncol(tempim),nrow=nrow(tempim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(data.stamp[,1])-(x.pix[i]-imsxl),y=1:length(data.stamp[1,])-(y.pix[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x.pix-x.pix[i]+1),y=(y.pix-y.pix[i]+1), pch=3)
        #Convert pix values onto asinh scale
        #X-Axis Limits; focus on 0 with ± limits around there
        #Stretch data to a useful scale
        stretchscale=ifelse(lim<0,1.5*10^(abs(lim)+1),1.5*10^(-1*(lim-1)))
        #Apply transformation
        tempvecstretch<-magmap(tempvec,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        #Calculate Axes Major Tick Mark Labels
        axespoints<-c(-10^(lim:(lim-2)),0,10^((lim-2):lim))
        #Calculate Axes Major Tick Mark Lenghts
        asinhtcls<-c(rep(0.5,length(axespoints)),rep(0.2,(length(axespoints)-1)*4))
        #Add Minor tick marks
        axespoints<-c(axespoints,c(-10^(lim:(lim-2)),10^((lim-2):lim))/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*2/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*3/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*4/5)
        #Calculate Tick locations on transformed axes
        asinhticks=magmap(axespoints,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',stretchscale=stretchscale)$map
        #Remove overlapping Tick Marks
        if (any(duplicated(asinhticks))) {
          while(length(which(duplicated(asinhticks)))>0) {
            kk<-which(duplicated(asinhticks))[1]
            indall<-which(asinhticks==asinhticks[kk])
            ind2<-indall[which.min(abs(axespoints[indall]))]
            axespoints<-c(axespoints[ind2],axespoints[-1*indall])
            asinhticks<-c(asinhticks[ind2],asinhticks[-1*indall])
            asinhtcls<-c(asinhtcls[ind2],asinhtcls[-1*indall])
          }
        }
        #Get pixel histogram on transformed axes
        pix<-hist(as.numeric(tempvecstretch),plot=FALSE,breaks=seq(0,1,length=100))
        #Plot Histogram and Count axis
        magplot(x=rev(rev(pix$breaks)[-1]),y=pix$counts,type='s',xlab='Randoms Pixel Values (pix)',ylab="Count",side=2,xlim=c(0,1),ylim=c(1,10^(log10(max(pix$counts))+1)),main="Pixel Histogram",log='y')
        #Convert Labels to Pretty style
        labs<-floor(log10(abs(axespoints)))
        pref<-ifelse(axespoints<0,"-","")
        labs<-paste(pref,ifelse(labs>0,"1e+","1e"),labs,"",sep="")
        labs[which(labs=="1e+Inf")]<-"0"
        labs[which(labs=="1e-Inf")]<-"0"
        check = grep("1e+", labs, fixed = TRUE)
        labs[check] = paste(sub("1e+", "10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("1e-", labs, fixed = TRUE)
        labs[check] = paste(sub("1e-", "10^{-", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e+", labs, fixed = TRUE)
        labs[check] = paste(sub("e+", "*x*10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e-", labs, fixed = TRUE)
        labs[check] = paste(sub("e-", "*x*10^{-", labs[check], fixed = TRUE), "}", sep = "")
        axespoints<-parse(text=labs)
        #Draw X major and minor ticks
        axis(1,asinhticks[which(asinhtcls==0.5)],labels=FALSE,tcl=0.5)
        axis(1,asinhticks[which(asinhtcls==0.2)],labels=FALSE,tcl=0.2)
        mtext('Blanks Pixel Values (pix)', 1, line = 2)
        #Draw major tick labels
        ind<-which(asinhtcls==0.5 & (labs=="0" | abs(asinhticks-0.5)>0.1))
        axis(1,asinhticks[ind],labels=axespoints[ind],tcl=0)
        #Draw histograms for each bin.
        for(iter in 1:numIters) {
          if (length(which(!is.na(tempvec[iter,])))>0) {
            tempvecstretch<-magmap(tempvec[iter,],lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
            tmp<-hist(tempvecstretch,plot=FALSE,breaks=pix$breaks)
            lines(x=rev(rev(pix$breaks)[-1]),y=tmp$counts,type='s',col=hsv(seq(2/3,0,length=numIters))[iter])
          }
        }
        abline(v=magmap(dat$randMean.mean/(max(sumap,na.rm=T)),lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Blanks Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",signif(dat$randMean.mean/max(sumap,na.rm=TRUE),digits=3)," (per pix)\nStd Dev = ",signif(dat$randMean.SD/sqrt(max(sumap,na.rm=TRUE)),digits=3)," (per pix)",sep=""),cex=0.6)
        boxx<-magmap(flux/sumap,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.2),boxwex=2)
        if (toFile) { dev.off() }
      }else{
        ranaps<-tempim*0
        tempvec<-matrix(NA,ncol=length(which(ap!=0)),nrow=numIters)
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          ranaps[as.matrix(ind)]<-ranaps[as.matrix(ind)]+gdapv*0.1
          tempvec[iter,]<-tempim[as.matrix(ind)]*gdapv
          #Sum shifted image and Aperture to return Flux
          val<-which(!is.na(tempim[as.matrix(ind)]))
          temp<-tempvec[iter,val]
          val<-gdapv[val]
          sumap[iter]<-sum(val)
          #Sum Flux, correcting for masked pixels
          flux[iter]<-sum(temp*val,na.rm=TRUE)
        }
        #Return Result
        s<-which(is.finite(flux)&is.finite(sumap))
        if (length(s)!=length(flux)) {
          flux<-flux[s]
          sumap<-sumap[s]
        }
        wflux<-sum(flux*sumap)/sum(sumap)
        wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
        wmad<-weightedMad(flux,sumap)
        if (sigclip>0) {
          for (iter in 1:nclip) {
            ind<-which((flux-wflux)/wsd <= sigclip)
            wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
            wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
            wmad<-weightedMad(flux[ind],sumap[ind])
          }
        }
        dat=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
        lim<-(ceiling(log10(max(abs(quantile(tempim,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        if (toFile) { CairoPNG(file=file.path(path,paste(cat.id[i],"_rancor.png",sep="")),height=6*res,width=10*res,res=res) }
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(data.stamp[,1])-(x.pix[i]-imsxl),y=1:length(data.stamp[1,])-(y.pix[i]-imsyl),matrix(magmap(tempim,lo=0.01,hi=0.99,stretch='asinh')$map,ncol=ncol(tempim),nrow=nrow(tempim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(data.stamp[,1])-(x.pix[i]-imsxl),y=1:length(data.stamp[1,])-(y.pix[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x.pix-x.pix[i]+1),y=(y.pix-y.pix[i]+1), pch=3)
        #Convert pix values onto asinh scale
        #X-Axis Limits; focus on 0 with ± limits around there
        #Stretch data to a useful scale
        stretchscale=ifelse(lim<0,1.5*10^(abs(lim)+1),1.5*10^(-1*(lim-1)))
        #Apply transformation
        tempvecstretch<-magmap(tempvec,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        #Calculate Axes Major Tick Mark Labels
        axespoints<-c(-10^(lim:(lim-2)),0,10^((lim-2):lim))
        #Calculate Axes Major Tick Mark Lenghts
        asinhtcls<-c(rep(0.5,length(axespoints)),rep(0.2,(length(axespoints)-1)*4))
        #Add Minor tick marks
        axespoints<-c(axespoints,c(-10^(lim:(lim-2)),10^((lim-2):lim))/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*2/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*3/5,c(-10^(lim:(lim-2)),10^((lim-2):lim))*4/5)
        #Calculate Tick locations on transformed axes
        asinhticks=magmap(axespoints,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',stretchscale=stretchscale)$map
        #Remove overlapping Tick Marks
        if (any(duplicated(asinhticks))) {
          while(length(which(duplicated(asinhticks)))>0) {
            kk<-which(duplicated(asinhticks))[1]
            indall<-which(asinhticks==asinhticks[kk])
            ind2<-indall[which.min(abs(axespoints[indall]))]
            axespoints<-c(axespoints[ind2],axespoints[-1*indall])
            asinhticks<-c(asinhticks[ind2],asinhticks[-1*indall])
            asinhtcls<-c(asinhtcls[ind2],asinhtcls[-1*indall])
          }
        }
        #Get pixel histogram on transformed axes
        pix<-hist(as.numeric(tempvecstretch),plot=FALSE,breaks=seq(0,1,length=100))
        #Plot Histogram and Count axis
        magplot(x=rev(rev(pix$breaks)[-1]),y=pix$counts,type='s',xlab='Randoms Pixel Values (pix)',ylab="Count",side=2,xlim=c(0,1),main="Pixel Histogram",log='y',ylim=c(1,10^(log10(max(pix$counts))+1)))
        #Convert Labels to Pretty style
        labs<-floor(log10(abs(axespoints)))
        pref<-ifelse(axespoints<0,"-","")
        labs<-paste(pref,ifelse(labs>0,"1e+","1e"),labs,"",sep="")
        labs[which(labs=="1e+Inf")]<-"0"
        labs[which(labs=="1e-Inf")]<-"0"
        check = grep("1e+", labs, fixed = TRUE)
        labs[check] = paste(sub("1e+", "10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("1e-", labs, fixed = TRUE)
        labs[check] = paste(sub("1e-", "10^{-", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e+", labs, fixed = TRUE)
        labs[check] = paste(sub("e+", "*x*10^{", labs[check], fixed = TRUE), "}", sep = "")
        check = grep("e-", labs, fixed = TRUE)
        labs[check] = paste(sub("e-", "*x*10^{-", labs[check], fixed = TRUE), "}", sep = "")
        axespoints<-parse(text=labs)
        #Draw X major and minor ticks
        axis(1,asinhticks[which(asinhtcls==0.5)],labels=FALSE,tcl=0.5)
        axis(1,asinhticks[which(asinhtcls==0.2)],labels=FALSE,tcl=0.2)
        mtext('Randoms Pixel Values (pix)', 1, line = 2)
        #Draw major tick labels
        ind<-which(asinhtcls==0.5 & (labs=="0" | abs(asinhticks-0.5)>0.1))
        axis(1,asinhticks[ind],labels=axespoints[ind],tcl=0)
        #Draw histograms for each bin.
        for(iter in 1:numIters) {
          if (length(which(!is.na(tempvec[iter,])))>0) {
            tempvecstretch<-magmap(tempvec[iter,],lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
            tmp<-hist(tempvecstretch,plot=FALSE,breaks=pix$breaks)
            lines(x=rev(rev(pix$breaks)[-1]),y=tmp$counts,type='s',col=hsv(seq(2/3,0,length=numIters))[iter])
          }
        }
        abline(v=magmap(dat$randMean.mean/(max(sumap,na.rm=T)),lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Random Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",signif(dat$randMean.mean/max(sumap,na.rm=TRUE),digits=3)," (per pix)\nStd Dev = ",signif(dat$randMean.SD/sqrt(max(sumap,na.rm=TRUE)),digits=3)," (per pix)",sep=""),cex=0.6)
        boxx<-magmap(flux/sumap,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.2),boxwex=2)
        if (toFile) { dev.off() }
      }
    }
  }
  return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
}

plotRanCor<-function(id_g,x_p,y_p,im_mask,imm_mask,ap_mask,stamplims,imstamplims,masklims,remask=TRUE,numIters=1E3,path="./",plotall=FALSE,sigclip=3,nclip=3,res=120){
  if(length(ap_mask) != length(stamplims[,1])){stop('ap_mask and stamplim lengths do not match!')}
  if (remask) { if(length(ap_mask) != length(masklims[,1])){stop('ap_mask and masklim lengths do not match!')} }
  cutup<-TRUE
  if(length(ap_mask) != length(im_mask)){
    if (is.matrix(im_mask)) {
      cutup<-FALSE
    } else {
      stop("ap_mask and im_mask lengths do not match!")
    }
  }
  if (remask) {
    if(length(ap_mask) != length(imm_mask)) {
      if (cutup){
        stop("ap_mask and imm_mask lengths do not match!")
      } else if (!is.matrix(imm_mask)) {
        stop("im_mask and imm_mask are not of the same form (both matrix or both list)")
      }
    }
  }
  if (!remask) {
    dir.create(file.path(path,"RanCorIms"),showWarnings=FALSE)
    path=file.path(path,"RanCorIms")
  } else {
    dir.create(file.path(path,"BlanksCorIms"),showWarnings=FALSE)
    path=file.path(path,"BlanksCorIms")
  }

  if(plotall) {
    rand<-1:length(x_p)
  } else {
    rand<-sample(length(x_p),min(10,length(x_p)))
  }
  if (cutup) {
    for (i in rand) {
     origim=im_mask[[i]]
    maskim=imm_mask[[i]]
         ap=ap_mask[[i]]
     sxl=stamplims[i,1]
     sxh=stamplims[i,2]
     syl=stamplims[i,3]
     syh=stamplims[i,4]
     imsxl=imstamplims[i,1]
     imsxh=imstamplims[i,2]
     imsyl=imstamplims[i,3]
     imsyh=imstamplims[i,4]

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
      if(remask){
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
        lim<-(floor(log10(max(abs(quantile(origim,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        CairoPNG(file=file.path(path,paste(id_g[i],"_blankscor.png",sep="")),height=6*res,width=10*res,res=res)
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(origim[,1])-(x_p[i]-imsxl),y=1:length(origim[1,])-(y_p[i]-imsyl),matrix(magmap(tempim,stretch='asinh')$map,ncol=ncol(tempim),nrow=nrow(tempim)),col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(origim[,1])-(x_p[i]-imsxl),y=1:length(origim[1,])-(y_p[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x_p-x_p[i]+1),y=(y_p-y_p[i]+1), pch=3)
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
        abline(v=magmap(dat$randMean.mean,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Blanks Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",round(dat$randMean.mean,digits=3),"\nStd Dev = ",round(dat$randMean.SD,digits=3),sep=""),cex=0.6)
        boxx<-magmap(flux,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.2),boxwex=2)
        dev.off()
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
        CairoPNG(file=file.path(path,paste(id_g[i],"_rancor.png",sep="")),height=6*res,width=10*res,res=res)
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(origim[,1])-(x_p[i]-imsxl),y=1:length(origim[1,])-(y_p[i]-imsyl),magmap(origim,stretch='asinh')$map,col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(origim[,1])-(x_p[i]-imsxl),y=1:length(origim[1,])-(y_p[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x_p-x_p[i]+1),y=(y_p-y_p[i]+1), pch=3)
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
        abline(v=magmap(dat$randMean.mean,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Random Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",round(dat$randMean.mean,digit=3),"\nStd Dev = ",round(dat$randMean.SD,digit=3),sep=""),cex=0.6)
        boxx<-magmap(flux,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))+0.5),boxwex=2)
        dev.off()
      }
    }
  } else {
    for (i in rand) {
         ap=ap_mask[[i]]
     sxl=stamplims[i,1]
     sxh=stamplims[i,2]
     syl=stamplims[i,3]
     syh=stamplims[i,4]
     imsxl=imstamplims[i,1]
     imsxh=imstamplims[i,2]
     imsyl=imstamplims[i,3]
     imsyh=imstamplims[i,4]

      #Get number of cols & rows in image stamp
      nc<-ncol(im_mask)
      nr<-nrow(im_mask)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      #Mask object aperture
      tempim<-im_mask
      tempim[sxl:sxh,syl:syh][which(zapsmall(ap)!=0,arr.ind=TRUE)]<-NA
      #Get Aperture details
      gdap<-which(ap!=0)
      gdapv<-ap[which(ap!=0)]
      #If Remasking
      if(remask){
        #Set mask 0s to NA
        imm_mask[which(imm_mask==0)]<-NA
        #Multiply Image and Mask together
        tempim<-tempim*imm_mask
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
        lim<-(floor(log10(max(abs(quantile(im_mask,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        CairoPNG(file=file.path(path,paste(id_g[i],"_blankscor.png",sep="")),height=6*res,width=10*res,res=res)
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(im_mask[,1])-(x_p[i]-imsxl),y=1:length(im_mask[1,])-(y_p[i]-imsyl),magmap(tempim,stretch='asinh')$map,col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(im_mask[,1])-(x_p[i]-imsxl),y=1:length(im_mask[1,])-(y_p[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x_p-x_p[i]+1),y=(y_p-y_p[i]+1), pch=3)
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
        abline(v=magmap(dat$randMean.mean,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Blanks Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",round(dat$randMean.mean,digit=3),"\nStd Dev = ",round(dat$randMean.SD,digit=3),sep=""),cex=0.6)
        boxx<-magmap(flux,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.2),boxwex=2)
        dev.off()
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
        lim<-(floor(log10(max(abs(quantile(tempim,c(0.001,0.999),na.rm=TRUE))))))
        if (!is.finite(lim)) { next }
        CairoPNG(file=file.path(path,paste(id_g[i],"_rancor.png",sep="")),height=6*res,width=10*res,res=res)
        layout(cbind(1,2))
        mar<-par("mar")
        par(mar=mar*c(1,0.8,1,0.2))
        image(x=1:length(im_mask[,1])-(x_p[i]-imsxl),y=1:length(im_mask[1,])-(y_p[i]-imsyl),magmap(tempim,stretch='asinh')$map,col=hsv(seq(2/3,0,length=256)),axes=FALSE,ylab="",xlab="",main="Image Stamp",asp=1,useRaster=TRUE)
        image(x=1:length(im_mask[,1])-(x_p[i]-imsxl),y=1:length(im_mask[1,])-(y_p[i]-imsyl),ranaps,col=hsv(0,0,0,alpha=0:100/100),add=TRUE,useRaster=TRUE)
        magaxis(side=1:4,labels=FALSE)
        magaxis(side=1:2,xlab="X (pix)",ylab="Y (pix)")
        points(x=(x_p-x_p[i]+1),y=(y_p-y_p[i]+1), pch=3)
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
        abline(v=magmap(dat$randMean.mean,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col=hsv(0,0,0,alpha=0.7),lty=1)
        abline(v=magmap(0,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map,col='darkgreen')
        legend('topright',legend=c("Random Flux; Mean"),col=hsv(0,0,0),lty=c(1),cex=0.6)
        label('topleft',lab=paste("Histograms show:\nBlack - All Randoms Pix\nColoured - Individual Randoms\nMean Est = ",round(dat$randMean.mean,digit=3),"\nStd Dev = ",round(dat$randMean.SD,digit=3),sep=""),cex=0.6)
        boxx<-magmap(flux,lo=-1*10^(lim),hi=10^(lim),range=c(0,1),type='num',stretch='asinh',clip='NA',stretchscale=stretchscale)$map
        boxplot(boxx,horizontal=TRUE,axes=FALSE,add=TRUE,pch=8,at=10^(log10(max(pix$counts))-1.2),boxwex=2)
        dev.off()
      }
    }
  }
  return=NULL
}

ran.cor<-function(data.stamp,ap.stamp,mask.stamp=NULL,ap.stamp.lims=NULL,data.stamp.lims=NULL,rem.mask=FALSE,numIters=1E2,mpi.opts="",sigclip=3,nclip=0,cat.x=NULL,cat.y=NULL,rand.x=NULL,rand.y=NULL,ran.main.mask.lim=0.99){
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

  if (is.null(cat.x)) {
    x.pix<-rowMeans(rbind(data.stamp.lims[,cbind(1,2)]))
  } else {
    if (length(cat.x) != length(data.stamp.lims[,1])) {
      stop("Supplied cat.x is not of the same length as the ap.stamp")
    }
    x.pix<-floor(cat.x)
  }
  if (is.null(cat.y)) {
    y.pix<-rowMeans(rbind(data.stamp.lims[,cbind(3,4)]))
  } else {
    if (length(cat.y) != length(data.stamp.lims[,1])) {
      stop("Supplied cat.y is not of the same length as the ap.stamp")
    }
    y.pix<-floor(cat.y)
  }

  if (cutup) {
    output<-foreach(origim=data.stamp, maskim=mask.stamp, ap=ap.stamp,
                    sxl=ap.stamp.lims[,1],sxh=ap.stamp.lims[,2],syl=ap.stamp.lims[,3],syh=ap.stamp.lims[,4],
                    catx=x.pix,caty=y.pix,.export=c('rand.x','rand.y'),.options.mpi=mpi.opts, .combine='rbind') %dopar% {
      #Get number of cols & rows in image stamp
      nc<-ncol(origim)
      nr<-nrow(origim)
      #Get shift numbers
      if (is.null(rand.x)) {
        dx<-round(runif(numIters, min=-1*floor(nr/2),max=floor(nr/2)))
      } else {
        if (length(rand.x != numIters)) {
          warning("Overwriting requested numIters with supplied number of rand.x entries")
          numIters<-length(rand.x)
        }
        dx<-catx-floor(rand.x)
      }
      if (is.null(rand.y)) {
        dy<-round(runif(numIters, min=-1*floor(nc/2),max=floor(nc/2)))
      } else {
        if (length(rand.y != numIters)) {
          warning("Overwriting requested numIters with supplied number of rand.y entries")
          numIters<-length(rand.y)
        }
        dy<-caty-floor(rand.y)
      }
      ind.use<-which(abs(dx) <= ceiling(nc/2) & abs(dy) <= ceiling(nr/2))
      if (length(ind.use) != numIters) {
        warning("There are supplied Randoms positions that are beyond the limits of the image\nThese are discarded.")
        dx<-dx[ind.use]
        dy<-dy[ind.use]
        numIters<-length(ind.use)
      }
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      sdap<-rep(NA,numIters)
      madap<-rep(NA,numIters)
      #If doing randoms, still Mask object aperture
      if (!rem.mask) {
        origim[sxl:sxh,syl:syh][which(zapsmall(ap)>=1-ran.main.mask.lim,arr.ind=TRUE)]<-NA
      }
      #Get aperture details
      gdap<-which(ap!=0)
      gdapv<-ap[which(ap!=0)]
      #If Remasking
      if(rem.mask){
        #Set mask 0s to NA
        maskim[which(maskim==0)]<-NA
        #Multiply Image and Mask together
        tempim<-origim*maskim
        #For Niters, calculate random Flux
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          temp<-tempim[as.matrix(ind)]*gdapv
          #Sum shifted image and Aperture to return Flux
          val<-which(!is.na(tempim[as.matrix(ind)]))
          temp<-temp[val]
          val<-gdapv[val]
          sumap[iter]<-sum(val)
          flux[iter]<-sum(temp*val,na.rm=TRUE)
        }
        #Return Result
        s<-which(is.finite(flux)&is.finite(sumap)&(sumap > 0))
        if (length(s)!=length(flux)) {
          flux<-flux[s]
          sumap<-sumap[s]
        }
        wflux<-sum(flux*sumap)/sum(sumap)
        wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
        wmad<-weightedMad(flux,sumap)
        if (nclip>0) {
          for (iter in 1:nclip) {
            ind<-which((flux-wflux)/wsd <= sigclip)
            wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
            wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
            wmad<-weightedMad(flux[ind],sumap[ind])
          }
        }
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
      }else{
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          temp<-origim[as.matrix(ind)]*gdapv
          #Sum shifted image and Aperture to return Flux
          val<-which(!is.na(origim[as.matrix(ind)]))
          temp<-temp[val]
          val<-gdapv[val]
          sumap[iter]<-sum(val)
          flux[iter]<-sum(temp*val,na.rm=TRUE)
        }
        #Return Result
        s<-which(is.finite(flux)&is.finite(sumap)&(sumap > 0))
        if (length(s)!=length(flux)) {
          flux<-flux[s]
          sumap<-sumap[s]
        }
        wflux<-sum(flux*sumap)/sum(sumap)
        wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
        wmad<-weightedMad(flux,sumap)
        if (nclip>0) {
          for (iter in 1:nclip) {
            ind<-which((flux-wflux)/wsd <= sigclip)
            wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
            wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
            wmad<-weightedMad(flux[ind],sumap[ind])
          }
        }
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux) & sumap>0)),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
      }
    }
  } else {
    #Make Temp Image
    if (rem.mask) {
      tempim<-mask.stamp
      #Set mask 0s to NA
      tempim[which(tempim==0)]<-NA
      #Multiply Image
      tempim<-tempim*data.stamp
    } else {
      tempim<-data.stamp
    }
    output<-foreach(ap=ap.stamp, .noexport=ls(envir=environment()), .export=c('rand.x','rand.y','tempim','numIters','rem.mask','sigclip','nclip'),
                    sxl=ap.stamp.lims[,1],sxh=ap.stamp.lims[,2],syl=ap.stamp.lims[,3],syh=ap.stamp.lims[,4],catx=x.pix,caty=y.pix,
                    .options.mpi=mpi.opts, .combine='rbind') %dopar% {
      #Get number of cols & rows in image stamp
      nc<-ncol(tempim)
      nr<-nrow(tempim)
      #Get shift numbers
      if (is.null(rand.x)) {
        dx<-round(runif(numIters, min=1,max=nr))
      } else {
        if (length(rand.x != numIters)) {
          warning("Overwriting requested numIters with supplied number of rand.x entries")
          numIters<-length(rand.x)
        }
        dx<-catx-floor(rand.x)
      }
      if (is.null(rand.y)) {
        dy<-round(runif(numIters, min=1,max=nc))
      } else {
        if (length(rand.y != numIters)) {
          warning("Overwriting requested numIters with supplied number of rand.y entries")
          numIters<-length(rand.y)
        }
        dy<-caty-floor(rand.y)
      }
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      sdap<-rep(NA,numIters)
      madap<-rep(NA,numIters)
      #Get good ap indicies
      gdap<-which(zapsmall(ap)!=0)
      gdapv<-ap[gdap]
      #For Niters, calculate random Flux
      for (iter in 1:numIters) {
        #Calculate Shift Indicies
        xind<-((1:(nr+1)+dx[iter])%%(nr+1))
        xind<-xind[xind>0][1:length(ap[,1])]
        yind<-((1:(nc+1)+dy[iter])%%(nc+1))
        yind<-yind[yind>0][1:length(ap[1,])]
        ind<-expand.grid(xind,yind)[gdap,]
        temp<-tempim[as.matrix(ind)]*gdapv
        #Sum shifted image and Aperture to return Flux
        val<-which(!is.na(tempim[as.matrix(ind)]))
        temp<-temp[val]
        val<-gdapv[val]
        sumap[iter]<-sum(val)
        #Sum Flux, correcting for masked pixels
        flux[iter]<-sum(temp*val,na.rm=TRUE)
      }
      #Return Result
      s<-which(is.finite(flux)&is.finite(sumap)&(sumap > 0))
      if (length(s)!=length(flux)) {
        flux<-flux[s]
        sumap<-sumap[s]
      }
      wflux<-sum(flux*sumap)/sum(sumap)
      wsd<-sqrt(sum(sumap * (flux - wflux)^2) * (sum(sumap)/(sum(sumap)^2 - sum(sumap))))
      wmad<-weightedMad(flux,sumap)
      if (nclip>0) {
        for (iter in 1:nclip) {
          ind<-which((flux-wflux)/wsd <= sigclip)
          wflux<-sum(flux[ind]*sumap[ind])/sum(sumap[ind])
          wsd<-sqrt(sum(sumap[ind] * (flux[ind] - wflux)^2) * (sum(sumap[ind])/(sum(sumap[ind])^2 - sum(sumap[ind]))))
          wmad<-weightedMad(flux[ind],sumap[ind])
        }
      }
      return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
    }
  }
  #Output the foreach data
  return=output
}

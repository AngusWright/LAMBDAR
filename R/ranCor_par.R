rancor.par<-function(im_mask,imm_mask,ap_mask,stamplims,masklims,remask=TRUE,numIters=1E3,mpiopts="",sigclip=3,nclip=0){
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

  if (cutup) {
    output<-foreach(origim=im_mask, maskim=imm_mask, ap=ap_mask,
                    sxl=stamplims[,1],sxh=stamplims[,2],syl=stamplims[,3],syh=stamplims[,4],
                    mxl=masklims[,1],mxh=masklims[,2],myl=masklims[,3],myh=masklims[,4], .options.mpi=mpiopts, .combine='rbind') %dopar% {
      #Get number of cols & rows in image stamp
      nc<-ncol(origim)
      nr<-nrow(origim)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      sdap<-rep(NA,numIters)
      madap<-rep(NA,numIters)
      #Mask object aperture
      origim[sxl:sxh,syl:syh][which(zapsmall(ap)!=0,arr.ind=TRUE)]<-NA
      #Get aperture details
      gdap<-which(ap!=0)
      gdapv<-ap[which(ap!=0)]
      #If Remasking
      if(remask){
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
        s<-which(is.finite(flux)&is.finite(sumap))
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
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
      }
    }
  } else {
    #Make Temp Image
    if (remask) {
      tempim<-imm_mask
      #Set mask 0s to NA
      tempim[which(tempim==0)]<-NA
      #Multiply Image
      tempim<-tempim*im_mask
    } else {
      tempim<-im_mask
    }
    output<-foreach(ap=ap_mask, .noexport=ls(envir=environment()), .export=c('tempim','numIters','remask'),
                    sxl=stamplims[,1],sxh=stamplims[,2],syl=stamplims[,3],syh=stamplims[,4],
                    mxl=masklims[,1],mxh=masklims[,2],myl=masklims[,3],myh=masklims[,4], .options.mpi=mpiopts, .combine='rbind') %dopar% {
      #Get number of cols & rows in image stamp
      nc<-ncol(tempim)
      nr<-nrow(tempim)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
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
      return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=TRUE),randAp.SD=sd(sumap,na.rm=TRUE),randAp.MAD=mad(sumap,na.rm=TRUE))
    }
  }
  #Output the foreach data
  return=output
}

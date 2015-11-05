rancor.par<-function(im_mask,imm_mask,ap_mask,stamplims,masklims,remask=TRUE,numIters=1E3,mpiopts=""){
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
        browser()
        stop("im_mask and imm_mask are not of the same form (both matrix or both list)")
      }
    }
  }

  if (cutup) {
    output<-foreach(origim=im_mask, maskim=imm_mask, ap=ap_mask,
                    sxl=stamplims[,1],sxh=stamplims[,2],syl=stamplims[,3],syh=stamplims[,4],
                    mxl=masklims[,1],mxh=masklims[,2],myl=masklims[,3],myh=masklims[,4], .options.mpi=mpiopts, .combine='rbind') %dopar% {
    #for (i in 1:length(ap_mask)) {
    # origim=im_mask[[i]]
    #maskim=imm_mask[[i]]
    #     ap=ap_mask[[i]]
    # sxl=stamplims[i,1]
    # sxh=stamplims[i,2]
    # syl=stamplims[i,3]
    # syh=stamplims[i,4]
    #  mxl=masklims[i,1]
    #  mxh=masklims[i,2]
    #  myl=masklims[i,3]
    #  myh=masklims[i,4]

      #Get number of cols & rows in image stamp
      nc<-ncol(origim)
      nr<-nrow(origim)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      #maskfrac<-rep(NA,numIters)
      #Mask object aperture
      origim[sxl:sxh,syl:syh][which(zapsmall(ap)!=0,arr.ind=T)]<-NA
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
          sumap[iter]<-sum(gdapv[which(!is.na(tempim[as.matrix(ind)]))])
          flux[iter]<-sum(temp,na.rm=T)/sumap[iter]
        }
        #Return Result
        wflux<-sum(flux*sumap,na.rm=T)/sum(sumap,na.rm=T)
        wsd<-sqrt(sum(sumap^2,na.rm=T)/sum(sumap,na.rm=T)^2*var(flux*sumap,na.rm=T))
        wmad<-sqrt(sum(sumap^2,na.rm=T)/(sum(sumap,na.rm=T)^2)*mad(flux*sumap,na.rm=T)^2)
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=T),randAp.SD=sd(sumap,na.rm=T),randAp.MAD=mad(sumap,na.rm=T))
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
          sumap[iter]<-sum(gdapv[which(!is.na(origim[as.matrix(ind)]))])
          flux[iter]<-sum(temp,na.rm=T)/sumap[iter]
        }
        #Return Result
        wflux<-sum(flux*sumap,na.rm=T)/sum(sumap,na.rm=T)
        wsd<-sqrt(sum(sumap^2,na.rm=T)/sum(sumap,na.rm=T)^2*var(flux*sumap,na.rm=T))
        wmad<-sqrt(sum(sumap^2,na.rm=T)/(sum(sumap,na.rm=T)^2)*mad(flux*sumap,na.rm=T)^2)
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=T),randAp.SD=sd(sumap,na.rm=T),randAp.MAD=mad(sumap,na.rm=T))
      }
    }
  } else {
    output<-foreach(ap=ap_mask, .noexport=ls(envir=environment()), .export=c('im_mask','imm_mask','numIters','remask'),
                    sxl=stamplims[,1],sxh=stamplims[,2],syl=stamplims[,3],syh=stamplims[,4],
                    mxl=masklims[,1],mxh=masklims[,2],myl=masklims[,3],myh=masklims[,4], .options.mpi=mpiopts, .combine='rbind') %dopar% {
    #for (i in 1:length(ap_mask)) {
    # origim=im_mask[[i]]
    #maskim=imm_mask[[i]]
    #     ap=ap_mask[[i]]
    # sxl=stamplims[i,1]
    # sxh=stamplims[i,2]
    # syl=stamplims[i,3]
    # syh=stamplims[i,4]
    #  mxl=masklims[i,1]
    #  mxh=masklims[i,2]
    #  myl=masklims[i,3]
    #  myh=masklims[i,4]

      #Get number of cols & rows in image stamp
      nc<-ncol(im_mask)
      nr<-nrow(im_mask)
      #Get shift numbers
      dx<-round(runif(numIters, min=1,max=nr))
      dy<-round(runif(numIters, min=1,max=nc))
      #Initialise Vectors
      flux<-rep(NA,numIters)
      sumap<-rep(NA,numIters)
      #Get good ap indicies
      gdap<-which(ap!=0)
      gdapv<-ap[which(ap!=0)]
      #If Remasking
      if(remask){
        #Set mask 0s to NA
        imm_mask[which(imm_mask==0)]<-NA
        #Multiply Image and Mask together
        tempim<-im_mask*imm_mask
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
          sumap[iter]<-sum(gdapv[which(!is.na(tempim[as.matrix(ind)]))])
          flux[iter]<-sum(temp,na.rm=T)/sumap[iter]
        }
        #Return Result
        wflux<-sum(flux*sumap,na.rm=T)/sum(sumap,na.rm=T)
        wsd<-sqrt(sum(sumap^2,na.rm=T)/sum(sumap,na.rm=T)^2*var(flux*sumap,na.rm=T))
        wmad<-sqrt(sum(sumap^2,na.rm=T)/(sum(sumap,na.rm=T)^2)*mad(flux*sumap,na.rm=T)^2)
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=T),randAp.SD=sd(sumap,na.rm=T),randAp.MAD=mad(sumap,na.rm=T))
      }else{
        for (iter in 1:numIters) {
          #Calculate Shift Indicies
          xind<-((1:(nr+1)+dx[iter])%%(nr+1))
          xind<-xind[xind>0][1:length(ap[,1])]
          yind<-((1:(nc+1)+dy[iter])%%(nc+1))
          yind<-yind[yind>0][1:length(ap[1,])]
          ind<-expand.grid(xind,yind)[gdap,]
          temp<-im_mask[as.matrix(ind)]*gdapv
          #Sum shifted image and Aperture to return Flux
          sumap[iter]<-sum(gdapv[which(!is.na(im_mask[as.matrix(ind)]))])
          flux[iter]<-sum(temp,na.rm=T)/sumap[iter]
        }
        #Return Result
        wflux<-sum(flux*sumap,na.rm=T)/sum(sumap,na.rm=T)
        wsd<-sqrt(sum(sumap^2,na.rm=T)/sum(sumap,na.rm=T)^2*var(flux*sumap,na.rm=T))
        wmad<-sqrt(sum(sumap^2,na.rm=T)/(sum(sumap,na.rm=T)^2)*mad(flux*sumap,na.rm=T)^2)
        return=data.frame(randMean.mean=wflux,randMean.SD=wsd,randMean.MAD=wmad,nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=T),randAp.SD=sd(sumap,na.rm=T),randAp.MAD=mad(sumap,na.rm=T))
      }
    }
  }
  #Output the foreach data
  return=output
}

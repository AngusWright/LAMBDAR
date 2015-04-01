rancor.par<-function(im_mask,imm_mask,ap_mask,stamplims,masklims,remask=TRUE,numIters=1E3,mpiopts=""){
  if(length(ap_mask) != length(im_mask)){stop('ap_mask and im_mask lengths do not much!')}
  if (remask) { if(length(ap_mask) != length(imm_mask)){stop('ap_mask and imm_mask lengths do not much!')} }
  if(length(ap_mask) != length(stamplims[,1])){stop('ap_mask and stamplim lengths do not much!')}
  if (remask) { if(length(ap_mask) != length(masklims[,1])){stop('ap_mask and masklim lengths do not much!')} }

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
    fluxmed<-rep(NA,numIters)
    #maskfrac<-rep(NA,numIters)
    #If Remasking
    if(remask){
      #Set mask 0s to NA
      maskim[which(maskim==0)]<-NA
      #Multiply Image and Mask together
      tempim<-origim*maskim
      #For Niters, calculate random Flux
      for (iter in 1:numIters) {
        #Calculate Shift Indicies
        #xind=(1:(nr+1)+dx[iter])%%(nr+1)
        #yind=(1:(nc+1)+dy[iter])%%(nc+1)
        temp<-tempim[(1:(nr+1)+dx[iter])%%(nr+1),(1:(nc+1)+dy[iter])%%(nc+1)][which(ap!=0)]
        #Sum shifted image and Aperture to return Flux
        flux[iter]<-mean(temp,na.rm=T)
        fluxmed[iter]<-median(temp,na.rm=T)
        #Sum shifted mask and Aperture, divided by sum Aperture, to return mask fraction
        #maskfrac[iter]<-1-sum(maskim[(1:(nr+1)+dx[iter])%%(nr+1),(1:(nc+1)+dy[iter])%%(nc+1)][which(ap!=0,arr.ind=T)])/length(which(ap!=0))
      }
      #Get number of unmasked randoms
      #noMask<-which(maskfrac==0)
      #Return Result
      #return=data.frame(randMean.mean=mean(flux),randMed.mean=mean(fluxmed),randMean.SD=sd(flux),randMed.SD=sd(fluxmed),randMean.noMask.mean=mean(flux[noMask]),randMean.noMask.SD=sd(flux[noMask]),nNoMask=length(noMask),meanMaskFrac=mean(maskfrac))
      return=data.frame(randMean.mean=mean(flux),randMed.mean=mean(fluxmed),randMean.SD=sd(flux),randMed.SD=sd(fluxmed))
    }else{
      for (iter in 1:numIters) {
        #Calculate Shift Indicies
        temp<-origim[(1:(nr+1)+dx[iter])%%(nr+1),(1:(nc+1)+dy[iter])%%(nc+1)][which(ap!=0)]
        #Sum shifted image and Aperture to return Flux
        flux[iter]<-mean(temp,na.rm=T)
        fluxmed[iter]<-median(temp,na.rm=T)
      }
      #Return Result
      #return=data.frame(randMean.mean=mean(flux),randMed.mean=mean(fluxmed),randMean.SD=sd(flux),randMed.SD=sd(fluxmed),randMean.noMask.mean=NA,randMean.noMask.SD=NA,nNoMask=NA,meanMaskFrac=NA)
      return=data.frame(randMean.mean=mean(flux),randMed.mean=mean(fluxmed),randMean.SD=sd(flux),randMed.SD=sd(fluxmed))
    }
  }
  #Output the foreach data
  return=output
}

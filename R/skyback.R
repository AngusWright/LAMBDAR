skyback<-function(ra,dec,cutlo=0,cuthi=100,origim,astrom,maskim,remmask=TRUE,radweight=1,clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,probcut=3){
  if(length(ra) != length(dec)){stop('ra and dec lengths do not much!')}
  if(length(ra) != length(cutlo)){stop('ra and cutlo lengths do not much!')}
  if(length(ra) != length(cuthi)){stop('ra and cuthi lengths do not much!')}
  cutlo[cutlo<hardlo*PSFFWHMinPIX]=hardlo*PSFFWHMinPIX
  cuthi[cuthi<hardhi*PSFFWHMinPIX]=hardhi*PSFFWHMinPIX
  ravec=ra;decvec=dec;cutlovec=round(cutlo);cuthivec=round(cuthi)
  pixlocvec=round(ad2xy(ra,dec,astrom))
  
  output=foreach(i=1:length(ra),.combine='rbind')%dopar%{
   #for (i in 1:length(ra)){
    ra=ravec[i];dec=decvec[i];cutlo=cutlovec[i];cuthi=cuthivec[i];pixloc=pixlocvec[i,]
    #Find extreme pixels to cut out
    xlocs=pixloc[1]+(-cuthi:cuthi)
    ylocs=pixloc[2]+(-cuthi:cuthi)
    #Find object location on new pixel grid
    xcen=cuthi+1-length(which(xlocs<0))
    ycen=cuthi+1-length(which(ylocs<0))
    #Select only pixels which are inside the image bounds
    xsel=which(xlocs>0 & xlocs<=astrom$NAXIS[1])
    ysel=which(ylocs>0 & ylocs<=astrom$NAXIS[2])
    #Trim to above
    xlocs=xlocs[xsel]
    ylocs=ylocs[ysel]
    #Create new cutout image, either the raw pixels, or multiplied through by the sourcemask
    if(remmask){
      tempim=origim$dat[[1]][xlocs,ylocs]*maskim$dat[[1]][xlocs,ylocs]
    }else{
      tempim=origim$dat[[1]][xlocs,ylocs]
    }
    #All ref pixels for new image
    tempref=as.matrix(expand.grid(1:length(xsel),1:length(ysel)))
    #Corresponding radii for new pixels from the object of interest
    temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
    #Keep only pixels inside the radius bounds given by cutlo and cuthi
    keep=temprad>cutlo & temprad<cuthi
    #Trim
    tempref=tempref[keep,]
    tempval=tempim[tempref]
    temprad=temprad[keep]
    #If sourcemask is used ignore pixels that exactly equal 0 (since these will belong to masked pixels)
    if(remmask){
      temprad=temprad[tempval!=0]
      tempval=tempval[tempval!=0]
    }
    #Do iterative <probcut>-sigma pixel clipping
    if(clipiters>0){
      probcut<-1-pnorm(probcut)
      for(i in 1:clipiters){
        vallims=2*median(tempval)-quantile(tempval,probcut)
        temprad=temprad[tempval<vallims]
        tempval=tempval[tempval<vallims]
      }
    }
    #Find the running medians for the data
    tempmedian=magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T)
    tempylims=tempmedian$ysd
    tempy=tempmedian$y
    #Calculate worst case sky error- the sd of the medians calculated
    skyerr=sd(tempy)
    #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
    weights=1/((tempmedian$x^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
    #Generate Sky RMS 
    skyRMS=as.numeric((quantile(tempval,0.5)-quantile(tempval,pnorm(-1))))
    #Determine Pearsons Test for Normality p-value for sky
    skyRMSpval=pearson.test(tempval)$p.value
    #Find the weighted mean of the medians
    sky=sum(tempy*weights)/(sum(weights))
    #Now we iterate until no running medians are outside the 1-sigma bound of the sky
    while(any(!(tempylims[,1]<=sky & tempylims[,2]>=sky)) & all(!(tempylims[,1]<=sky & tempylims[,2]>=sky))==FALSE){
      tempy=tempy[tempylims[,1]<=sky & tempylims[,2]>=sky]
      weights=weights[tempylims[,1]<=sky & tempylims[,2]>=sky]
      tempylims=rbind(tempylims[tempylims[,1]<=sky & tempylims[,2]>=sky,])
      sky=sum(tempy*weights)/(sum(weights))
    }
    #Find the number of running medians that agree with the final sky within error bounds (max=10)
    Nnearsky=length(which(tempylims[,1]<=sky & tempylims[,2]>=sky))
    #Organise data into data.frame for foreach
    data.frame(sky=sky,skyerr=skyerr,Nnearsky=Nnearsky,skyRMS=skyRMS,skyRMSpval=skyRMSpval)
  }
  #Output the foreach data
  return=output
}

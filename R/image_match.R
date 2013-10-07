image_match <-
function(ra0=-999,dec0=-999,pathroot="",inpim=NA,cutrad=1/60,outputfits=NA, match=NA){
    # Procedure takes a fits image and produces a new  
    # image, cropped to a region <cutrad> in diameter 
    # and centered around ra0 dec0
    # All inputs can be vectorised to length 'n', such that the output
    # will be a series of 'n' images whose names are the n 
    # dimensions of <fitsoutname>
    # NB:
    #   > ra0, dec0,fitsoutname must be of the same dimension  
    #   > If ra0, dec0, and fitsoutname have dimension > 0, and
    #     if cutrad, pathroot, inpim are of dimension 0, then 
    #     all ra0,dec0,fitsoutname are given the same cutrad, 
    #     pathroot, and inpim
    #   > If not supplied, cutrad == 1deg, pathroot==""
    #   > If not supplied, ra0 and dec0 will be set to 
    #   > the Image Centre
    

    # Ensure ra0, dec0, and fitsoutname have same dimensions
    len<-max(length(ra0),length(inpim))
    if((((length(dec0       )-len)>0)&(length(dec0    )!=1))| # Must be exactly the same length or 1
       (((length(ra0        )-len)>0)&(length(ra0     )!=1))| # Must be exactly the same length or 1
       (((length(cutrad     )-len)>0)&(length(cutrad  )!=1))| # Must be exactly the same length or 1
       (((length(pathroot   )-len)>0)&(length(pathroot)!=1))| # Must be exactly the same length or 1
       (((length(outputfits )-len)>0)&(all(!is.na(outputfits)))) | # Must be exactly the same length or absent
       (((length(inpim      )-len)>0)&(length(inpim   )!=1))) # Must be exactly the same length or 1
    { stop("Mismatched Vector Dimensions in Inputs") }
     
    match<-unique(match)

    # Check Other input dimensions
    if (len>1) {
      if (length(ra0     )==1) { ra0     <-replicate(len, ra0     ) }
      if (length(dec0    )==1) { dec0    <-replicate(len, dec0    ) }
      if (length(pathroot)==1) { pathroot<-replicate(len, pathroot) }
      if (length(cutrad  )==1) { cutrad  <-replicate(len, cutrad  ) }
      if (length(inpim   )==1) { inpim   <-replicate(len, inpim   ) }
    }  

    # File definitions
    image <- paste(pathroot, inpim, sep="")
        
    # Setup Arrays
    xcen<-array(NA, len)
    ycen<-array(NA, len)
    naxis1<-array(NA, len)
    naxis2<-array(NA, len)
    pixscale1<-array(NA, len)
    pixscale2<-array(NA, len)
    size<-array(NA, len)
    type<-array(NA, len)
    signed<-array(NA, len)
    headbytes<-array(NA, len)
    # Get header values
    for (i in 1:len) {
      astr_struc<-read.astr(image[i])
      naxis1[i]<-astr_struc$NAXIS[1]
      naxis2[i]<-astr_struc$NAXIS[2]
      pixscale1[i]<-astr_struc$CD[1,1]
      pixscale2[i]<-astr_struc$CD[2,2]
      # If needed, set ra0 and dec0
      if (ra0[i] ==-999) { ra0 <-replicate(len,as.numeric(astr_struc$CRVAL[1])) }
      if (dec0[i]==-999) { dec0<-replicate(len,as.numeric(astr_struc$CRVAL[2])) }
      # Calculate pixel positions
      pixpos <- ad2xy(ra0[i],dec0[i],astr_struc)
      xcen[i] <- pixpos[,"X"]
      ycen[i] <- pixpos[,"Y"]
      # Get data precision type
      switch(astr_struc$BITPIX, "-64" = {
          size[i] = 8           # 64-bit float (double)
          type[i] = 'numeric'
          signed[i]=TRUE
      }, "-32" = {
          size[i] = 4           # 32-bit float (single)
          type[i] = 'numeric'
          signed[i]=TRUE
      }, "32" = {
          size[i] = 4           # 32-bit signed int (double)
          type[i] = 'integer'
          signed[i]=TRUE
      }, "16" = {
          size[i] = 2           # 16-bit signed int (single)
          type[i] = 'integer'
          signed[i]=TRUE
      }, "8" = {
          size[i] = 1           # 8-bit unsigned int
          type[i] = 'integer'
          signed[i]=FALSE
      }, stop("Unknown BITPIX request in .read.fits.image"))
      # Check that cutrad is sensible
      if (cutrad[i] < 0 ) { 
        warning("Cut Radius is < 0. Using 1deg radius instead")
        cutrad <- 1 
      }
    }

    # Define lower and upper limits (make sure the limits 
    # are within the image)
    cutradxpix <- ceiling(abs((cutrad/pixscale1))/2)*2+1
    cutradypix <- ceiling(abs((cutrad/pixscale2))/2)*2+1

    interpolate=FALSE
    if ((any(!is.na(match)))&(any(c(abs(pixscale1[match]),abs(pixscale2[match])) != abs(pixscale1[match][1])))) { 
      interpolate=TRUE 
      lores<-match[which(abs(pixscale1[match])!=min(abs(pixscale1[match])))]
      #print(paste("Lores Index:", lores))
      #lores<-which(abs(pixscale1)!=min(abs(pixscale1)))[1]
      cutradxpix[lores]<-cutradxpix[lores]+4
      cutradypix[lores]<-cutradypix[lores]+4
    }

    #Make Limits
    lowerx <- floor(xcen)-cutradxpix
    upperx <- floor(xcen)+cutradxpix
    lowery <- floor(ycen)-cutradypix
    uppery <- floor(ycen)+cutradypix
    cutxcen <- xcen-lowerx
    cutycen <- ycen-lowery
    # Adjust reference pixel location in cases where image
    # width changes
    cutxcen[which(lowerx<1)] <- abs(cutrad[which(lowerx<1)]/pixscale1[which(lowerx<1)])-(1-lowerx[which(lowerx<1)])
    cutycen[which(lowery<1)] <- abs(cutrad[which(lowery<1)]/pixscale2[which(lowery<1)])-(1-lowery[which(lowery<1)])
    # Change image width if pixel indexes are -ve or > naxis1/2
    lowerx[which(lowerx<1)]<-1
    #print(lowery)
    lowery[which(lowery<1)]<-1
    #print(lowery)
    upperx[which(upperx>naxis1)]<-naxis1[which(upperx>naxis1)]
    uppery[which(uppery>naxis2)]<-naxis2[which(uppery>naxis2)]

    #Check for Full Images
    fullim<-((lowerx==1)&(upperx==naxis1)&(lowery==1)&(uppery==naxis2))


    ncol<-upperx-lowerx+1
    nrow<-uppery-lowery+1
    #print(upperx)
    #print(lowerx)
    #print(uppery)
    #print(lowery)
    #print(ncol)
    #print(nrow)
    #print(cutrad)
    #print(cutradxpix)
    #print(cutradypix)
    #print(ncol*pixscale1/2)
    #print(nrow*pixscale2/2)
    
    #if ((any(ncol != ncol[1]))|(any(nrow != nrow[1]))) { stop("Input Images have mis-matched limits") }

    imlist<-list()
    skypos<-list()
    for (i in 1:len) {
      #print(paste("Image",i)) 
      # Output the subimage header to file
      header<-read.fitshdr(image[i])
      pipe<-file(image[i], 'rb')
      h<-make.fits.header(header)
      headbytes<-nchar(h, type='bytes')
      # Make/check skypos array
      skypos <- c(skypos, list(xy2radec(x=seq(1.5,ncol[i]+0.5,by=1),y=seq(1.5,nrow[i]+0.5,by=1),
                                        ra0=ra0[i],dec0=dec0[i],x0=cutxcen[i],y0=cutycen[i],
                                        xscale=pixscale1[i],yscale=pixscale2[i])))
      #if (any(!is.finite(skypos[[i]]))){print("NA produced in ra calc")}
      #print(paste(min(skypos[[i]][,'RA']),max(skypos[[i]][,'RA'])))
      #print(paste(min(skypos[[i]][,'DEC']),max(skypos[[i]][,'DEC'])))
      imlist$x[[i]]<-cbind(RA=skypos[[i]][,'RA']) 
      imlist$y[[i]]<-cbind(DEC=skypos[[i]][,'DEC'])
      # Use seek to get to the start of the subimage
      #if (((grepl("fuv",image[i],ignore.case=TRUE))|
           #(grepl("nuv",image[i],ignore.case=TRUE)))&
           #(!grepl("g15",image[i],ignore.case=TRUE))) { imstrt=naxis1[i]*((lowery[i]))+lowerx[i]+1 }
      #else { imstrt=naxis1[i]*((lowery[i]+1))+lowerx[i] }
      imstrt=naxis1[i]*((lowery[i]+1))+lowerx[i]
      seek(pipe, where=headbytes+(imstrt)*size[i],origin='start')
      if (!fullim[i]) {
        # Read the subimage with readBin, read each row, then go to the start 
        # of the next row
        subim<-array(0,dim=c(nrow[i],ncol[i]))
        for (j in 1:nrow[i]){
          subim[j,]<-readBin(pipe,type[i],n=ncol[i],size[i],signed=signed[i],endian='big')
          seek(pipe, where=(naxis1[i]-ncol[i])*size[i], origin='current')
        }
        subim<-(subim-min(subim))
        subim<-subim/max(subim)
        imlist$dat<-c(imlist$dat,list(subim))
        close(pipe)
      } else {
        imlist$dat<-c(imlist$dat,list(read.fits(image[i],hdu=0)$dat[[1]]))
      }
      if (all(!is.na(outputfits))) {
        header<-read.fitshdr(inpim[i])
        # New header vals
        header[which(header[,'key']=='NAXIS1'),'value']<-paste(ncol[i])
        header[which(header[,'key']=='NAXIS2'),'value']<-paste(nrow[i])
        header[which(header[,'key']=='CRVAL1'),'value']<-paste(ra0[i])
        header[which(header[,'key']=='CRVAL2'),'value']<-paste(dec0[i])
        header[which(header[,'key']=='CRPIX1'),'value']<-paste(cutxcen[i])
        header[which(header[,'key']=='CRPIX2'),'value']<-paste(cutycen[i])
        #header[which(header[,'key']=='CD1_1'),'value']<-as.numeric(header[which(header[,'key']=='CD1_1'),'value'])*-1
        #header[which(header[,'key']=='CD2_2'),'value']<-as.numeric(header[which(header[,'key']=='CD2_2'),'value'])*-1
        out<-matrix(imlist$dat[[i]], ncol=ncol[i], byrow=TRUE)
	write.fits(list(hdr=list(header),dat=list(out)), file=outputfits[i])
      }
    }
    #print("Done")
    if (interpolate) {
      #print("Interpolating")
      hires<-match[which(abs(pixscale1[match])==min(abs(pixscale1[match])))][1]
      #print(paste("HiRes is", hires))
      #hires<-which(abs(pixscale1)==min(abs(pixscale1)))[1]
      #print("Expanding grid")
      #print(length(skypos[[hires]][,'RA']))
      #print(length(skypos[[hires]][,'DEC']))
      grid<-expand.grid(skypos[[hires]][,'RA'],skypos[[hires]][,'DEC'])
      #plot(grid, xlim=c(min(grid[,1])-diff(range(grid[,1]))/2,max(grid[,1])+diff(range(grid[,1]))/2),ylim=c(min(grid[,2])-diff(range(grid[,2]))/2,max(grid[,2])+diff(range(grid[,2]))/2))
      for (i in lores) {
         #print("Assigning grid")
         imlist$x[[i]]<-cbind(RA=skypos[[hires]][,'RA'])
         imlist$y[[i]]<-cbind(DEC=skypos[[hires]][,'DEC'])
	#print(paste("Interpolating image",i))
        #points(expand.grid(skypos[[i]][,'RA'],skypos[[i]][,'DEC']), col=i, pch=i+1)
        imlist$dat[[i]]<-matrix(interp2D(x=grid[,1],y=grid[,2],list(x=skypos[[i]][,'RA'],y=skypos[[i]][,'DEC'],z=imlist$dat[[i]])),ncol=ncol[hires])
      }
      #points(ra0, dec0, col='red', lw=2.0, pch="#")
    }

    return=imlist

}


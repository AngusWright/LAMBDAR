crop_im <-
function(ra0=-999,dec0=-999,pathroot="",inpim=NA,cutrad=1,fitsoutname=NA){
#Details {{{
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
#   > the Image Centre }}}

    cat(paste("   Cropping Image",inpim,"   "))
    # Ensure ra0, dec0, and fitsoutname have same dimensions {{{
    if (ra0!=-999) { len<-length(ra0) } else { len<-length(inpim) } # Length is determined by ra0 vec or inpim vec
    if((((length(dec0       )-len)!=0)&(length(dec0    )!=1))| # Must be exactly the same length or 1
       (((length(ra0        )-len)!=0)&(length(ra0     )!=1))| # Must be exactly the same length or 1
        ((length(fitsoutname)-len)!=0)|                        # Must be exactly the same length
       (((length(cutrad     )-len)!=0)&(length(cutrad  )!=1))| # Must be exactly the same length or 1
       (((length(pathroot   )-len)!=0)&(length(pathroot)!=1))| # Must be exactly the same length or 1
       (((length(inpim      )-len)!=0)&(length(inpim   )!=1))) # Must be exactly the same length or 1
    { stop("Mismatched Vector Dimensions in Inputs") }
    #}}}

    # Check Other input dimensions {{{
    if (len>1) {
      if (length(ra0     )==1) { ra0     <-replicate(len, ra0     ) }
      if (length(dec0    )==1) { dec0    <-replicate(len, dec0    ) }
      if (length(pathroot)==1) { pathroot<-replicate(len, pathroot) }
      if (length(cutrad  )==1) { cutrad  <-replicate(len, cutrad  ) }
      if (length(inpim   )==1) { inpim   <-replicate(len, inpim   ) }
    }#}}}

    # File definition {{{
    image <- paste(pathroot, inpim, sep="")
    #}}}

    # Setup Arrays {{{
    xcen<-array(NA, len)
    ycen<-array(NA, len)
    naxis1<-array(NA, len)
    naxis2<-array(NA, len)
    pixsize<-array(NA, len)
    size<-array(NA, len)
    type<-array(NA, len)
    signed<-array(NA, len)
    headbytes<-array(NA, len)
    #}}}

    # Get header & cut radius values {{{
    for (i in 1:len) {
      #Get astrometry values {{{
      astr_struc<-read.astr(image[i])
      naxis1[i]<-astr_struc$NAXIS[1]
      naxis2[i]<-astr_struc$NAXIS[2]
      pixsize[i]<-astr_struc$CD[1,1]
      #}}}
      # If needed, set ra0 and dec0 {{{
      if (ra0[i] ==-999) { ra0[i] <-as.numeric(astr_struc$CRVAL[1]) }
      if (dec0[i]==-999) { dec0[i]<-as.numeric(astr_struc$CRVAL[2]) }
      #}}}
      # Calculate pixel positions {{{
      pixpos <- ad2xy(ra0[i],dec0[i],astr_struc)
      xcen[i] <- pixpos[,"X"]
      ycen[i] <- pixpos[,"Y"]
      #}}}
      # Get data precision type {{{
      switch(astr_struc$BITPIX, "-64" = {
          #64-bit float (double) {{{
          size[i] = 8
          type[i] = 'numeric'
          signed[i]=TRUE
          #}}}
      }, "-32" = {
          #32-bit float (single) {{{
          size[i] = 4
          type[i] = 'numeric'
          signed[i]=TRUE
          #}}}
      }, "32" = {
          #32-bit signed int (double) {{{
          size[i] = 4
          type[i] = 'integer'
          signed[i]=TRUE
          #}}}
      }, "16" = {
          #16-bit signed int (single) {{{
          size[i] = 2
          type[i] = 'integer'
          signed[i]=TRUE
          #}}}
      }, "8" = {
          #8-bit unsigned int {{{
          size[i] = 1
          type[i] = 'integer'
          signed[i]=FALSE
          #}}}
      }, stop("Unknown BITPIX request in .read.fits.image"))
      #}}}
      # Check that cutrad is sensible {{{
      if (cutrad[i] < 0 ) {
        warning("Cut Radius is < 0. Using 1deg radius instead")
        cutrad <- 1
      }#}}}
    } #}}}

    # Define lower and upper limits {{{
    cutradpix <- abs((cutrad/pixsize))
    lowerx <- floor(xcen-cutradpix)
    upperx <- floor(xcen+cutradpix)
    lowery <- floor(ycen-cutradpix)
    uppery <- floor(ycen+cutradpix)
    cutxcen <- abs(cutrad/pixsize)+1.5
    cutycen <- abs(cutrad/pixsize)+1.5
    #}}}

    #If pixel limits are beyond image limits... {{{
    #Shift the cut-centre {{{
    cutxcen[which(lowerx<1)] <- abs(cutrad[which(lowerx<1)]/pixsize[which(lowerx<1)])+1.5-(1-lowerx[which(lowerx<1)])
    cutycen[which(lowery<1)] <- abs(cutrad[which(lowery<1)]/pixsize[which(lowery<1)])+1.5-(1-lowery[which(lowery<1)])
    #}}}
    #Change the limits {{{
    lowerx[which(lowerx<1)]<-1
    lowery[which(lowery<1)]<-1
    upperx[which(upperx>naxis1)]<-naxis1[which(upperx>naxis1)]
    uppery[which(uppery>naxis2)]<-naxis2[which(uppery>naxis2)]
    #}}}
    #}}}

    #Define Row & Col lengths {{{
    ncol<-upperx-lowerx+1
    nrow<-uppery-lowery+1
    #}}}

    #Crop Image(s) {{{
    for (i in 1:len) {
      #Get the Image header {{{
      header<-read.fitshdr(inpim[i])
      #}}}
      #Update header values for new image {{{
      header[which(header[,'key']=='NAXIS1'),'value']<-paste(ncol[i])
      header[which(header[,'key']=='NAXIS2'),'value']<-paste(nrow[i])
      header[which(header[,'key']=='CRPIX1'),'value']<-paste(as.numeric(header[which(header[,'key']=='CRPIX1'),'value'])-lowerx)
      header[which(header[,'key']=='CRPIX2'),'value']<-paste(as.numeric(header[which(header[,'key']=='CRPIX2'),'value'])-lowery)
      #}}}
      #Write the new header to file {{{
      pipe<-file(inpim[i], 'rb')
      pipeout<-file(fitsoutname[i], 'wb')
      h<-make.fits.header(header)
      writeChar(h, pipeout, eos=NULL)
      headbytes<-nchar(h, type='bytes')
      #}}}
      # Use seek to get to the start of the subimage {{{
      imstrt=naxis1[i]*(lowery[i]-1)+lowerx[i]-1
      seek(pipe, where=headbytes+(imstrt)*size[i],origin='start')
      #}}}
      #Read the subimage, row by row, and immediately output to new image {{{
      for (j in 1:nrow[i]){
        writeBin(readBin(pipe,type[i],n=ncol[i],size[i],signed=signed[i],endian='big'),pipeout,size=size[i],endian='big')
        seek(pipe, where=(naxis1[i]-ncol[i])*size[i], origin='current')
      }
      #}}}
      #Close Pipes {{{
      close(pipe)
      close(pipeout)
      #}}}
    }#}}}
    cat(" - Done\n")
    return=NULL
}


crop.fits.image <-
function(ra0=-999,dec0=-999,path.root="./",inpim=NA,crop.radius=1,fitsoutname=NA,data.extn=1){
#Details {{{
# Procedure takes a fits image and produces a new
# image, cropped to a region <crop.radius> in diameter
# and centered around ra0 dec0
# All inputs can be vectorised to length 'n', such that the output
# will be a series of 'n' images whose names are the n
# dimensions of <fitsoutname>
# NB:
#   > ra0, dec0,fitsoutname must be of the same dimension
#   > If ra0, dec0, and fitsoutname have dimension > 0, and
#     if crop.radius, path.root, inpim are of dimension 0, then
#     all ra0,dec0,fitsoutname are given the same crop.radius,
#     path.root, and inpim
#   > If not supplied, crop.radius == 1deg, path.root==""
#   > If not supplied, ra0 and dec0 will be set to
#   > the Image Centre }}}

    cat(paste("   Cropping Image",inpim,"   "))
    # Ensure ra0, dec0, and fitsoutname have same dimensions {{{
    if (ra0!=-999) { len<-length(ra0) } else { len<-length(inpim) } # Length is determined by ra0 vec or inpim vec
    if((((length(dec0       )-len)!=0)&(length(dec0    )!=1))| # Must be exactly the same length or 1
       (((length(ra0        )-len)!=0)&(length(ra0     )!=1))| # Must be exactly the same length or 1
        ((length(fitsoutname)-len)!=0)|                        # Must be exactly the same length
       (((length(crop.radius     )-len)!=0)&(length(crop.radius  )!=1))| # Must be exactly the same length or 1
       (((length(path.root   )-len)!=0)&(length(path.root)!=1))| # Must be exactly the same length or 1
       (((length(inpim      )-len)!=0)&(length(inpim   )!=1))) # Must be exactly the same length or 1
    { sink(type='message'); stop("Mismatched Vector Dimensions in Inputs") }
    #}}}

    # Check Other input dimensions {{{
    if (len>1) {
      if (length(ra0     )==1) { ra0     <-replicate(len, ra0     ) }
      if (length(dec0    )==1) { dec0    <-replicate(len, dec0    ) }
      if (length(path.root)==1) { path.root<-replicate(len, path.root) }
      if (length(crop.radius  )==1) { crop.radius  <-replicate(len, crop.radius  ) }
      if (length(inpim   )==1) { inpim   <-replicate(len, inpim   ) }
    }#}}}

    # File definition {{{
    image <- paste(path.root, inpim, sep="")
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
      astr.struc<-read.astrometry(image[i],hdu=data.extn)
      if (all(is.na(astr.struc$CD))) { sink(type='message'); stop("FITS extension does not have a valid WCS astrometry") }
      naxis1[i]<-astr.struc$NAXIS[1]
      naxis2[i]<-astr.struc$NAXIS[2]
      pixsize[i]<-astr.struc$CD[1,1]
      #}}}
      # If needed, set ra0 and dec0 {{{
      if ((ra0[i]==-999)||(dec0[i]==-999)) { locs<-xy.to.ad(floor(naxis1[i]/2), floor(naxis2[i]/2), astr.struc) }
      if (ra0[i] ==-999) { ra0[i] <-locs[1] }
      if (dec0[i]==-999) { dec0[i]<-locs[2] }
      if (ra0[i] ==-999) { ra0[i] <-as.numeric(astr.struc$CRVAL[1]) }
      if (dec0[i]==-999) { dec0[i]<-as.numeric(astr.struc$CRVAL[2]) }
      #}}}
      # Calculate pixel positions {{{
      pixpos <- ad.to.xy(ra0[i],dec0[i],astr.struc)
      xcen[i] <- pixpos[,"X"]
      ycen[i] <- pixpos[,"Y"]
      #}}}
      # Get data precision type {{{
      switch(astr.struc$BITPIX, "-64" = {
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
      # Check that crop.radius is sensible {{{
      if (crop.radius[i] < 0 ) {
        warning("Cut Radius is < 0. Using 1deg radius instead")
        crop.radius <- 1
      }#}}}
    } #}}}

    # Define lower and upper limits {{{
    crop.radiuspix <- abs((crop.radius/pixsize))
    lowerx <- floor(xcen-crop.radiuspix)
    upperx <- floor(xcen+crop.radiuspix)
    lowery <- floor(ycen-crop.radiuspix)
    uppery <- floor(ycen+crop.radiuspix)
    #}}}

    #If pixel limits are beyond image limits... {{{
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
      header<-read.fitshdr(file.path(path.root,inpim[i]),hdu=data.extn)
      #}}}
      #Update header values for new image {{{
      header[which(header[,'key']=='NAXIS1'),'value']<-paste(ncol[i])
      header[which(header[,'key']=='NAXIS2'),'value']<-paste(nrow[i])
      header[which(header[,'key']=='CRPIX1'),'value']<-paste(as.numeric(header[which(header[,'key']=='CRPIX1'),'value'])+(1-lowerx[i]))
      header[which(header[,'key']=='CRPIX2'),'value']<-paste(as.numeric(header[which(header[,'key']=='CRPIX2'),'value'])+(1-lowery[i]))
      #}}}
      #Write the new header to file {{{
      pipe<-file(file.path(path.root,inpim[i]), 'rb')
      pipeout<-file(file.path(path.root,fitsoutname[i]), 'wb')
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
      #Pad with NULL to end on a 2880 boundary {{{
      nbytes<-ncol[i]*nrow[i]*size[i]+headbytes
      npad<-ceiling(nbytes/2880)*2880 - nbytes
      writeBin(raw(npad),pipeout,endian='big')
      #}}}
      #Close Pipes {{{
      close(pipe)
      close(pipeout)
      #}}}
    }#}}}
    cat(" - Done\n")
    return=NULL
}


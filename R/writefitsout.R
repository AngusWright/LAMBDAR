writefitsout <-
function(filename,image,header_in,xsize=NA,ysize=NA,xcent=NA,ycent=NA,pixeloffset=NA, nochange=FALSE,verbose=FALSE,diagnostic=FALSE){
  #Write out a fits file with associated astrometry
  #xcent and ycent are the position in the input image of what will become the central pixel in the output image
  #these are in the IDL convention where the first pixel is (0,0)
  #nochange means there is no change to the astrometry in the fits header

  #Initialise headers
  header<-header_in
  astr_out<-NULL

  if (!(nochange)) {
    if ((!(is.finite(xsize)))|(!(is.finite(ysize)))|(!(is.finite(xcent)))|(!(is.finite(ycent)))) {
      sink(type="message")
      stop(paste('Incorrect Calling Syntax. Please Use:\n',
                 ' if nochange==TRUE:\n',
                 '                  > writefitsout(filename,image,header_in,<pixeloffset>, nochange=TRUE)\n',
                 ' if nochange==FALSE:\n',
                 '                  > writefitsout(filename,image,header_in,xsize,ysize,xcent,ycent,<pixeloffset>, <nochange==FALSE>)\n',
                 '                                                                            NB: commands inside <...> are optional'))
    }

    #Create new FITS header astrometry
    if (!(is.na(header["PC1_1" ,]))) { astr_out<-c(astr_out, list(PC1_1 =header["PC1_1" ,])) }
    if (!(is.na(header["PC1_2" ,]))) { astr_out<-c(astr_out, list(PC1_2 =header["PC1_2" ,])) }
    if (!(is.na(header["CDELT1",]))) { astr_out<-c(astr_out, list(CDELT1=header["CDELT1",])) }
    if (!(is.na(header["CDELT2",]))) { astr_out<-c(astr_out, list(CDELT2=header["CDELT2",])) }
    if (!(is.na(header["CRPIX1",]))) { astr_out<-c(astr_out, list(CRPIX1=header["CRPIX1",])) }
    if (!(is.na(header["CRPIX2",]))) { astr_out<-c(astr_out, list(CRPIX2=header["CRPIX2",])) }
    if (!(is.na(header["CRVAL1",]))) { astr_out<-c(astr_out, list(CRVAL1=header["CRVAL1",])) }
    if (!(is.na(header["CRVAL2",]))) { astr_out<-c(astr_out, list(CRVAL2=header["CRVAL2",])) }
    if (!(is.na(header["CD1_1" ,]))) { astr_out<-c(astr_out, list(CD1_1 =header["CD1_1" ,])) }
    if (!(is.na(header["CD2_1" ,]))) { astr_out<-c(astr_out, list(CD2_1 =header["CD2_1" ,])) }
    if (!(is.na(header["CD1_2" ,]))) { astr_out<-c(astr_out, list(CD1_2 =header["CD1_2" ,])) }
    if (!(is.na(header["CD2_2" ,]))) { astr_out<-c(astr_out, list(CD2_2 =header["CD2_2" ,])) }
    if (!(is.na(header["CROTA1",]))) { astr_out<-c(astr_out, list(CROTA1=header["CROTA1",])) }
    if (!(is.na(header["CROTA2",]))) { astr_out<-c(astr_out, list(CROTA2=header["CROTA2",])) }

    if ((is.finite(xcent)) & (is.finite(ycent))) {

      radec<-xy2ad(xcent,ycent,header)
      ra_deg<-radec[,"RA"]
      dec_deg<-radec[,"DEC"]


      #keep astr same as in astr_in except for the following tags:
      #      $NAXIS - 2 element array giving image size
      astr_out$NAXIS1<-xsize
      astr_out$NAXIS2<-ysize
      #      $CRPIX - 2 element double vector giving X and Y coordinates of reference
      #               pixel (def = NAXIS/2) in FITS convention (first pixel is 1,1)
      astr_out$CRPIX1<-1+floor(xsize/2.)
      astr_out$CRPIX2<-1+floor(ysize/2.)
      #      $CRVAL - 2 element double precision vector giving R.A. and DEC of
      #             reference pixel in DEGREES
      astr_out$CRVAL1<-ra_deg
      astr_out$CRVAL2<-dec_deg

    } else if (is.finite(pixeloffset)) {
      astr_out$CRPIX1<-astr_out$CRPIX1+pixeloffset
      astr_out$CRPIX2<-astr_out$CRPIX2+pixeloffset
    }

    if (!(is.null(astr_out$PC1_1 ))) { header["PC1_1" ,]<-astr_out$PC1_1 }
    if (!(is.null(astr_out$PC1_2 ))) { header["PC1_2" ,]<-astr_out$PC1_2 }
    if (!(is.null(astr_out$CDELT1))) { header["CDELT1",]<-astr_out$CDELT1}
    if (!(is.null(astr_out$CDELT2))) { header["CDELT2",]<-astr_out$CDELT2}
    if (!(is.null(astr_out$CRPIX1))) { header["CRPIX1",]<-astr_out$CRPIX1}
    if (!(is.null(astr_out$CRPIX2))) { header["CRPIX2",]<-astr_out$CRPIX2}
    if (!(is.null(astr_out$CRVAL1))) { header["CRVAL1",]<-astr_out$CRVAL1}
    if (!(is.null(astr_out$CRVAL2))) { header["CRVAL2",]<-astr_out$CRVAL2}
    if (!(is.null(astr_out$CD1_1 ))) { header["CD1_1" ,]<-astr_out$CD1_1 }
    if (!(is.null(astr_out$CD2_1 ))) { header["CD2_1" ,]<-astr_out$CD2_1 }
    if (!(is.null(astr_out$CD1_2 ))) { header["CD1_2" ,]<-astr_out$CD1_2 }
    if (!(is.null(astr_out$CD2_2 ))) { header["CD2_2" ,]<-astr_out$CD2_2 }
    if (!(is.null(astr_out$CROTA1))) { header["CROTA1",]<-astr_out$CROTA1}
    if (!(is.null(astr_out$CROTA2))) { header["CROTA2",]<-astr_out$CROTA2}

  }

  #-----Diagnostic-----#
  if (verbose) { message(paste("Writing FITS to file:",filename)) }
  if (diagnostic) {
    message(paste("Header to be output:"))
    sink(file=sinkfile, type="output", append=TRUE)
    print(cbind(rownames(header),header[,1]))
    sink(file=NULL, type="output")
  }
  #Output new FITS image with header
  headout<-cbind(rownames(header),header[,1])
  colnames(headout)<-c("key","value")
  fits<-list(hdr=list(headout), dat=list(image))
  write.fits(fits,file=filename,type="double",hdu=0)
  return=NULL
}

readpsf <-
function (env=NULL, filename,asperpix,apsize,stampmult,normalize=FALSE,gauss_fwhm_as=0, outenv=NULL) {

  message('--------------------------Read_PSF-------------------------------------')
 
  # Load Parameter Space
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env } 
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))
  

  if (gauss_fwhm_as > 0) {
    if (verbose) { message("Creating Gaussian PSF from input parameters") }
    psffwhm<-gauss_fwhm_as
    stampsizepix=(floor((ceiling(apsize/asperpix)+ceiling(psffwhm*stampmult/asperpix))/2)*2+1)
    psfsizepix=(((ceiling(psffwhm*stampmult/asperpix))/2)*2+1)

    sigma_as=gauss_fwhm_as/(2.*sqrt(2.*log(2.)))
    sigma_p=sigma_as/asperpix
    x0=floor(psfsizepix/2.)
    y0=floor(psfsizepix/2.)
    if (verbose) { message(paste("gauss_fwhm_as, x0, y0, height, sigma_p, sigma_as, psfsizepix, stampsizepix\n",round(gauss_fwhm_as, digits=3), x0, y0, height, sigma_p, sigma_as, psfsizepix, stampsizepix)) }
    
    zp<-foreach (xi=1:psfsizepix, .combine='cbind') %:% 
      foreach (yi=1:psfsizepix, .combine='c') %do% { 
        height*exp(-(((xi-x0)^2/(2*sigma_p^2))+((yi-y0)^2/(2*sigma_p^2))))
      }

    im_psf<-array(as.double(0), dim=c(stampsizepix,stampsizepix))
    im_psf[1:psfsizepix,1:psfsizepix]<-zp  # will be normalized to height

  } else {
    if (verbose) { message(paste("Reading PSF from file",filename)) }
    psf<-read.fits(filename,hdu=0)
    im_psf<-psf$dat[[1]]
    hdr_psf<-psf$hdr[[1]]
    if (normalize) { im_psf<-height*im_psf/max(im_psf) } 

    px<-(length(im_psf[,1]))
    py<-(length(im_psf[1,]))  
    #Get PSF FWHM
    xfwhm<-diff(range(which(im_psf[,py/2]>=max(im_psf)/2)))
    yfwhm<-diff(range(which(im_psf[px/2,]>=max(im_psf)/2)))
    fwhm_pix<-ceiling(mean(xfwhm,yfwhm))
    psffwhm<-fwhm_pix*asperpix
    stampsizepix<-(floor((ceiling(apsize/asperpix)+ceiling(fwhm_pix*stampmult))/2)*2+1)
    if (!(is.finite(stampsizepix)&is.finite(px)&is.finite(py))){ print(paste(px,py,xfwhm,yfwhm,fwhm_pix,asperpix,apsize,stampsizepix)) }
    if (min(px,py) < stampsizepix) {
      im_big<-array(0, dim=c(stampsizepix,stampsizepix))
      #NB: When using psfconvolve(), the psf location in the image 
      #    is ignored - therefore location of psf in im_big is irrelevant.
      ########################################################################
      # >>>IF NOT USING psfcolvolve(): LOCATION OF PSF IS HIGHLY RELEVANT<<< #
      ########################################################################
      im_big[1:px,1:py]<-im_psf
      im_psf<-im_big
    }
    if (verbose) { message(paste("image_fwhm_as, px, py, height, xfwhm_pix, yfwhm_pix, stampsizepix\n",round(psffwhm, digits=3), px, py, max(im_psf), xfwhm, yfwhm, stampsizepix)) }
    #Ensure PSF is centred on stamp by convolving with deltafunction
    delt<-array(0, dim=dim(im_psf))
    delt[ceiling(px/2),ceiling(py/2)]<-1
    oldsum<-sum(im_psf)
    im_psf<-convolvepsf(im_psf,delt,nomag2=TRUE)
    #Check nothing has changed except position
    if (sum(im_psf)!=oldsum) { message(paste("NOTE: PSF sum has changed from in realignment by",abs(oldsum-sum(im_psf))*100/oldsum,"%")) }
  }
  if (verbose) { message(paste("PSF FWHM =",round(psffwhm, digits=3),"arcsec")) 
                 message(paste("Stamp Size =",stampsizepix,"pixels")) }
  sumpsf<-as.single(sum(im_psf))
  #Parse Parameter Space
  detach(env)
  assign("psffwhm" , psffwhm  , envir = outenv)
  assign("sumpsf"  , sumpsf   , envir = outenv)
  
  im_psf
}

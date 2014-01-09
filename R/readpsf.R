readpsf <-
function (env=NULL, filename,asperpix,apsize,confidence,normalize=TRUE,gauss_fwhm_as=0, outenv=NULL) {

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
    psffwhm.pix<-gauss_fwhm_as/asperpix
    psfsigma.pix<-psffwhm.pix/(2*sqrt(2*log(2)))
    nsig<-sqrt(qchisq(confidence,df=2)) #Determine nsigma for desired confidence using chisq distribution
    psf.clip<-nsig*psfsigma.pix
    psf.clip<-psf.clip*2 # convert psf.clip from radius to diameter
    psfwidth<-psf.clip*asperpix
    stampsizepix=(floor((ceiling(apsize/asperpix)+ceiling(psf.clip/2))/2)*2+1)

    x0=floor(psf.clip/2.)
    y0=floor(psf.clip/2.)
    if (verbose) { message(paste("gauss_fwhm_as, x0, y0, sigma_p, stampsizepix\n",
                           round(gauss_fwhm_as, digits=3), x0, y0, psfsigma.pix, stampsizepix)) }

    zp<-foreach (xi=1:psf.clip, .combine='cbind') %:%
      foreach (yi=1:psf.clip, .combine='c') %do% {
        exp(-(((xi-x0)^2/(2*psfsigma.pix^2))+((yi-y0)^2/(2*psfsigma.pix^2))))
      }

    im_psf<-array(as.double(0), dim=c(stampsizepix,stampsizepix))
    im_psf[1:psf.clip,1:psf.clip]<-zp  # will be normalized to Unity

  } else {
    if (verbose) { message(paste("Reading PSF from file",filename)) }
    psf<-read.fits(filename,hdu=0)
    im_psf<-psf$dat[[1]]
    px<-length(im_psf[,1])
    py<-length(im_psf[1,])
    hdr_psf<-psf$hdr[[1]]
    #Get PSF FWHM
    psf.clip<-get.confidence(im_psf,confidence)
    psf.clip<-psf.clip*2 # convert psf.clip from radius to diameter
    psfwidth<-psf.clip*asperpix
    stampsizepix<-(floor((ceiling(apsize/asperpix)+ceiling(psf.clip/2))/2)*2+1)

    if (any(dim(im_psf) < stampsizepix)) {
      im_big<-array(0, dim=c(stampsizepix,stampsizepix))
      im_big[1:px,1:py]<-im_psf
      im_psf<-im_big
    }
    if (verbose) { message(paste("psf_width_as, psf.clip, height, stampsizepix\n",round(psfwidth, digits=3), psf.clip, max(im_psf), stampsizepix)) }
    #Ensure PSF is centred on stamp by convolving with deltafunction
    delt<-array(0, dim=dim(im_psf))
    delt[ceiling(px/2),ceiling(py/2)]<-1
    oldsum<-sum(im_psf)
    im_psf<-convolvepsf(im_psf,delt,nomag2=TRUE)
    #Check nothing has changed except position
    if (sum(im_psf)!=oldsum) { message(paste("NOTE: PSF sum has changed from in realignment by",abs(oldsum-sum(im_psf))*100/oldsum,"%")) }
  }
  if (verbose) { message(paste("PSF WIDTH =",round(psfwidth, digits=3),"arcsec"))
                 message(paste("Stamp Size =",stampsizepix,"pixels")) }
  if (normalize) {
    message("PSF is being normalised to Unity")
    im_psf<-im_psf/max(im_psf)
  }
  sumpsf<-as.single(sum(im_psf))
  #Parse Parameter Space
  detach(env)
  assign("psf.clip" , psf.clip  , envir = outenv)
  assign("psfwidth" , psfwidth  , envir = outenv)
  assign("sumpsf"   , sumpsf    , envir = outenv)

  im_psf
}

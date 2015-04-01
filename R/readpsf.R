readpsf <-
function (outenv=parent.env(environment()), filename,asperpix,apsize,confidence,normalize=TRUE,gauss_fwhm_as=0, env=NULL) {

  message('--------------------------Read_PSF-------------------------------------')

  # Load Parameter Space {{{
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  if(is.null(outenv)&!is.null(env)) { outenv<-env }
  else if (is.null(outenv)) {
    warning("Output Environment cannot be NULL; using parent env")
    outenv<-parent.env(environment())
  }
  #}}}

  #Read/Generate PSF {{{
  if (gauss_fwhm_as > 0) {
    #Generate Gaussian PSF {{{
    if (verbose) { message("Creating Gaussian PSF from input parameters") }
    #Get gaussian parameters {{{
    psffwhm.pix<-gauss_fwhm_as/asperpix
    psfsigma.pix<-psffwhm.pix/(2*sqrt(2*log(2)))
    if (confidence==1) {
      warning("Cannot have PSF Confidence = Unity when analytically deriving PSF; Using Maximum double value (1-1E-16)")
      confidence=1-1E-16
    }
    nsig<-sqrt(qchisq(confidence,df=2)) #Determine nsigma for desired confidence using chisq distribution
    psf.clip<-ceiling(nsig*psfsigma.pix)
    psf.clip<-psf.clip*2+1 # convert psf.clip from radius to diameter (and make sure it's odd)
    psfwidth<-psf.clip*asperpix
    stampsizepix=(floor((ceiling(defbuff*apsize*2/asperpix)+ceiling(psf.clip))/2)*2+5)

    x0=ceiling(psf.clip/2.)
    y0=ceiling(psf.clip/2.)
    if (verbose) { message(paste("gauss_fwhm_as, x0, y0, sigma_p, stampsizepix\n",
                           round(gauss_fwhm_as, digits=3), x0, y0, psfsigma.pix, stampsizepix)) }
    #}}}
    #Generate PSF pixel values {{{
    xi=matrix(1:psf.clip, nrow=psf.clip, ncol=psf.clip)
    yi=t(xi)
    zp<-exp(-(((xi-x0)^2/(2*psfsigma.pix^2))+((yi-y0)^2/(2*psfsigma.pix^2))))
    #}}}
    #Make PSF 2D array {{{
    im_psf<-array(as.double(0), dim=c(stampsizepix,stampsizepix))
    im_psf[1:psf.clip,1:psf.clip]<-zp  # will be normalized to Unity
    #}}}
    #}}}
  } else {
    #Read PSF from File {{{
    if (verbose) { message(paste("Reading PSF from file",filename)) }
    psf<-read.fits(filename,hdu=0)
    im_psf<-psf$dat[[1]]
    hdr_psf<-read.astr(filename)
    #}}}
    #Check that psf pixel scale is the same as the image {{{
    hdr_psf$CD<-abs(hdr_psf$CD)
    if (is.na(hdr_psf$CD[1,1])) {
      warning("No astrometry present in PSF Header. ASSUMING EQUIVALENT PIXEL SCALES")
      hdr_psf$CD[1:2,1:2]<-asperpix/3600
    }
    if ((hdr_psf$CD[1,1]*3600!=asperpix)||(hdr_psf$CD[2,2]*3600!=asperpix)) {
      #If it isn't, reinterpolate the psf onto the same pixel spacing
      narcsec_x<-((hdr_psf$NAXIS[1]-1)*hdr_psf$CD[1,1]*3600)
      narcsec_y<-((hdr_psf$NAXIS[2]-1)*hdr_psf$CD[2,2]*3600)
      x_oldres<-seq(0,narcsec_x,by=hdr_psf$CD[1,1]*3600)#+0.5*hdr_psf$CD[1,1]*3600
      y_oldres<-seq(0,narcsec_y,by=hdr_psf$CD[2,2]*3600)#+0.5*hdr_psf$CD[2,2]*3600
      x_newres<-seq(0,narcsec_x,by=asperpix            )#+0.5*asperpix
      y_newres<-seq(0,narcsec_y,by=asperpix            )#+0.5*asperpix
      grid<-expand.grid(x_newres,y_newres)
      im_psf<-matrix(interp2D(x=grid[,1],y=grid[,2],list(x=x_oldres,y=y_oldres,z=im_psf)),ncol=length(x_newres))
    }
    #}}}
    #Get PSF FWHM {{{
    px<-length(im_psf[,1])
    py<-length(im_psf[1,])
    psfvals<-rev(sort(im_psf))
    tempsum<-cumsum(psfvals)
    tempfunc<-approxfun(tempsum,psfvals)
    psfLimit<-tempfunc(confidence*max(tempsum, na.rm=TRUE))
    psf.clip<-diff(range(which(im_psf>=psfLimit, arr.ind=TRUE)[,1]))
    psf.clip<-max(psf.clip,diff(range(which(im_psf>=psfLimit, arr.ind=TRUE)[,2])))
    psfwidth<-psf.clip*asperpix
    stampsizepix<-(floor((ceiling(defbuff*apsize*2/asperpix)+ceiling(psf.clip))/2)*2+5)
    #}}}
    #If Needed, pad the PSF with 0s {{{
    if (any(dim(im_psf) < stampsizepix)) {
      im_big<-array(0, dim=c(stampsizepix,stampsizepix))
      im_big[1:px,1:py]<-im_psf
      im_psf<-im_big
    }
    #}}}
    if (verbose) { message(paste("psf_width_as, psf.clip, height, stampsizepix\n",round(psfwidth, digits=3), psf.clip, max(im_psf), stampsizepix)) }
    ##Ensure PSF is centred on a pixel centre by convolving with deltafunction {{{
    #delt<-array(0, dim=dim(im_psf))
    #delt[ceiling(px/2),ceiling(py/2)]<-1
    #oldsum<-sum(im_psf)
    #im_psf<-convolvepsf(im_psf,delt,mag2=FALSE)
    ##}}}
    ##Check nothing has changed except position {{{
    #if (sum(im_psf)!=oldsum) { message(paste("NOTE: PSF sum has changed from in realignment by",abs(oldsum-sum(im_psf))*100/oldsum,"%")) }
    ##}}}
  }
  #}}}
  #Notify {{{
  if (verbose) { message(paste("PSF WIDTH =",round(psfwidth, digits=3),"arcsec"))
                 message(paste("Stamp Size =",stampsizepix,"pixels")) }
  #}}}
  #Normalise {{{
  if (normalize) {
    message("PSF is being normalised to Unity")
    im_psf<-im_psf/max(im_psf)
  }
  #}}}
  #Finally, Remove any Negative Values in the PSF, and calculate the sum {{{
  im_psf[which(im_psf<0)]<-0
  sumpsf<-as.single(sum(im_psf))
  #}}}
  #Parse Parameter Space & Return{{{
  if (!is.null(env)) { detach(env) }
  assign("psf.clip" , psf.clip  , envir = outenv)
  assign("psfwidth" , psfwidth  , envir = outenv)
  assign("sumpsf"   , sumpsf    , envir = outenv)
  return=im_psf
  #}}}
}

read.psf <-
function (outenv=parent.env(environment()), filename,arcsec.per.pix,apsize,confidence,normalize=TRUE,gauss.fwhm.arcsec=0, env=NULL) {

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
  if (gauss.fwhm.arcsec != 0) {
    #Generate Gaussian PSF {{{
    if (gauss.fwhm.arcsec < 0) {
      message("Reading Gaussian PSF FWHM from datamap header")
      if (grepl('sig',psf.label.type,ignore.case=TRUE)) {
        #Add factor to convert pixel to arcsec, if needed
        if (grepl('pix',psf.label.type,ignore.case=TRUE)) { fact=arcsec.per.pix } else { fact=1 }
        gauss.sigma.arcsec<-as.numeric(read.fitskey(file=paste(path.root,path.work,data.map,sep=""),key=psf.label,hdu=ifelse(data.extn<=1,1,data.extn)))*fact
        gauss.fwhm.arcsec<-gauss.sigma.arcsec*(2*sqrt(2*log(2)))
        message(paste0("PSF width in header is as SIGMA: value = ",gauss.sigma.arcsec," arcsec (factor was: ",fact,")"))
      } else {
        #Add factor to convert pixel to arcsec, if needed
        if (grepl('pix',psf.label.type,ignore.case=TRUE)) { fact=arcsec.per.pix } else { fact=1 }
        gauss.fwhm.arcsec<-as.numeric(read.fitskey(file=paste(path.root,path.work,data.map,sep=""),key=psf.label,hdu=ifelse(data.extn<=1,1,data.extn)))*fact
        message(paste0("PSF width in header is as FWHM: value = ",gauss.fwhm.arcsec," arcsec (factor was: ",fact,")")) 
      }
      if (is.na(gauss.fwhm.arcsec)) { sink(type="message") ; stop("Read of PSF FWHM from datamap failed. Cannot continue") }
      if (verbose) { message(" -Done\n") }
    }
    if (verbose) { message("Creating Gaussian PSF from input parameters") }
    #Get gaussian parameters {{{
    psffwhm.pix<-gauss.fwhm.arcsec/arcsec.per.pix
    psfsigma.pix<-psffwhm.pix/(2*sqrt(2*log(2)))
    if (confidence==1) {
      warning("Cannot have PSF Confidence = Unity when analytically deriving PSF; Using Maximum double value (1-1E-16)")
      confidence=1-1E-16
    }
    nsig<-sqrt(qchisq(confidence,df=2)) #Determine nsigma for desired confidence using chisq distribution
    psf.clip<-ceiling(nsig*psfsigma.pix)
    psf.clip<-psf.clip*2+1 # convert psf.clip from radius to diameter (and make sure it's odd)
    psfwidth<-psf.clip*arcsec.per.pix
    if (psf.filt) {
      if (aperture.type==2) { 
        stampsizepix<-(floor(ceiling(def.buff*nsig*apsize/arcsec.per.pix)+(ceiling(psf.clip)/2))*2+5)
      } else {
        stampsizepix<-(floor(ceiling(def.buff*apsize/arcsec.per.pix)+(ceiling(psf.clip)/2))*2+5)
      }
    } else {
      if (aperture.type==2) { 
        stampsizepix<-(floor(ceiling(def.buff*nsig*apsize/arcsec.per.pix))*2+5)
      } else {
        stampsizepix<-(floor(ceiling(def.buff*apsize/arcsec.per.pix))*2+5)
      }
    }
    #if (psf.filt) {
    #  stampsizepix=(floor(ceiling(def.buff*apsize/arcsec.per.pix)+(ceiling(psf.clip)/2))*2+5)
    #} else {
    #  stampsizepix=(floor(ceiling(def.buff*apsize/arcsec.per.pix))*2+5)
    #}
    if (psf.clip > stampsizepix) {
      message(paste("WARNING: psf stamp is larger than the largest aperture stamp.\n",
      "It is being cut down in size, which will bias the aperture correction.\n",
      "To avoid this, increase the defBuff size in the parameter file."))
      psf.clip<-stampsizepix
      psfwidth<-psf.clip*arcsec.per.pix
    }
    x0=ceiling(psf.clip/2.)
    y0=ceiling(psf.clip/2.)
    if (verbose) { message(paste("gauss.fwhm.arcsec, x0, y0, sigma_p, stampsizepix\n",
                           round(gauss.fwhm.arcsec, digits=3), x0, y0, psfsigma.pix, stampsizepix)) }
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
    #Check psf file exists {{{
    if (!file.exists(filename)) {
      sink(type='message')
      stop("PSF Image does not exist at location specified:",filename)
    }#}}}
    psf<-try(read.fits(filename,hdu=0),silent=TRUE)
    if (class(psf)=="try-error") {
      sink(type='message')
      stop("PSF file read failed")
    }
    im_psf<-psf$dat[[1]]
    hdr_psf<-read.astrometry(filename)
    #}}}
    #Check that psf pixel scale is the same as the image {{{
    hdr_psf$CD<-abs(hdr_psf$CD)
    if (is.na(hdr_psf$CD[1,1])&is.na(hdr_psf$CD[2,2])) {
      warning("No astrometry present in PSF Header. ASSUMING EQUIVALENT PIXEL SCALES")
      hdr_psf$CD[1:2,1:2]<-arcsec.per.pix/3600
    } else if (is.na(hdr_psf$CD[1,1])) {
      warning("Incomplete astrometry present in PSF Header. Assuming Symmetry.")
      hdr_psf$CD[1,1]<-hdr_psf$CD[2,2]
    } else if (is.na(hdr_psf$CD[2,2])) {
      warning("Incomplete astrometry present in PSF Header. Assuming Symmetry.")
      hdr_psf$CD[2,2]<-hdr_psf$CD[1,1]
    }
    if ((abs(1-(hdr_psf$CD[1,1]*3600/arcsec.per.pix))>1E-3)||(abs(1-(hdr_psf$CD[2,2]*3600/arcsec.per.pix)>1E-3))) {
      if (hdr_psf$CD[1,1]*3600 > 10*arcsec.per.pix) { 
        #The psf is significantly lower resolution than the image; 
        #this is probably a mistake in the header, so ignore it
        message("WARNING: ignoring the PSF header information because it doesn't make sense...")
        message("         The PSF resolution is more than 10x lower than the image! ")
      } else { 
        #If it isn't, reinterpolate the psf onto the same pixel spacing
        narcsec_x<-((hdr_psf$NAXIS[1]-1)*hdr_psf$CD[1,1]*3600)
        narcsec_y<-((hdr_psf$NAXIS[2]-1)*hdr_psf$CD[2,2]*3600)
        x_oldres<-seq(0,narcsec_x,by=hdr_psf$CD[1,1]*3600)
        y_oldres<-seq(0,narcsec_y,by=hdr_psf$CD[2,2]*3600)
        x_newres<-seq(0,narcsec_x,by=arcsec.per.pix            )
        y_newres<-seq(0,narcsec_y,by=arcsec.per.pix            )
        grid<-expand.grid(x_newres,y_newres)
        im_psf<-matrix(interp.2d(x=grid[,1],y=grid[,2],list(x=x_oldres,y=y_oldres,z=im_psf))[,3],ncol=length(x_newres))
      }
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
    psfwidth<-psf.clip*arcsec.per.pix
    if (psf.filt) {
      stampsizepix=(floor(ceiling(def.buff*apsize/arcsec.per.pix)+(ceiling(psf.clip)/2))*2+5)
    } else {
      stampsizepix=(floor(ceiling(def.buff*apsize/arcsec.per.pix))*2+5)
    }
    #}}}
    #If Needed, pad the PSF with 0s {{{
    if (any(dim(im_psf) < stampsizepix)) {
      im_big<-array(0, dim=c(stampsizepix,stampsizepix))
      im_big[1:px,1:py]<-im_psf
      im_psf<-im_big
    }
    #}}}
    if (verbose) { message(paste("psf_width_as, psf.clip, height, stampsizepix\n",round(psfwidth, digits=3), psf.clip, max(im_psf), stampsizepix)) }
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
  psf.cen<-list(which(im_psf==max(im_psf),arr.ind=TRUE))
  #}}}
  #Finally, Remove any Negative Values in the PSF, and calculate the sum {{{
  im_psf[which(im_psf<0)]<-0
  sumpsf<-as.single(sum(im_psf))
  #}}}
  #Parse Parameter Space & Return{{{
  if (!is.null(env)) { detach(env) }
  assign("psf.cen"  , psf.cen   , envir = outenv)
  assign("psf.clip" , psf.clip  , envir = outenv)
  assign("psfwidth" , psfwidth  , envir = outenv)
  assign("sumpsf"   , sumpsf    , envir = outenv)
  return=im_psf
  #}}}
}

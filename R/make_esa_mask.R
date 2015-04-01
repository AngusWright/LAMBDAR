make_esa_mask <-
function(outenv=parent.env(environment()), env=NULL,ObsParm,padGals,col.corr=0,confuse=FALSE){
#Procedure creates exponential galaxy profiles, using input parameters
#from the catalogue, and places them (in order) onto stamps
#Procedure is parallelised to allow scaleability

  message("Creating Simulated Galaxies from Input Catalogue")
  message("Assumptions are:\nCatalogue Galaxy Apertures are 2.5*Kron Apertures\nCatalogue Fluxweights are Jy fluxes for use in Simulation")
  if (!quiet) { cat('Make_ESA_Mask   ') }
  message('--------------------------Make_ESA_Mask---------------------------------')

  # Load Parameter Space
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  if(is.null(outenv)&!is.null(env)) { outenv<-env }
  else if (is.null(outenv)) {
    warning("Output Environment cannot be NULL; using parent env")
    outenv<-parent.env(environment())
  }
  #

  #Convert semi-major (Kron) axis in arcsec to Major Effective Radius in pixels
  Reff_pix<-(a_g)*(1/(2.5*1.19))/(asperpix)
  #

  #Aperture Axis Ratio in [0,1]
  axrat<-b_g/a_g
  #

  #Get input Magnitudes
  inputmags<- -2.5*log10(fluxweight)+8.9
  #

  #Remove Galaxies with non-finite input Magnitudes
  index<-which(is.finite(inputmags))
  inputmags<-inputmags[index]
  a_g<-a_g[index]
  b_g<-b_g[index]
  ra_g<-ra_g[index]
  dec_g<-dec_g[index]
  x_g<-x_g[index]
  y_g<-y_g[index]
  theta_g<-theta_g[index]
  Reff_pix<-Reff_pix[index]
  fluxweight<-fluxweight[index]
  axrat<-axrat[index]
  id_g<-id_g[index]
  #

  if (padGals) {
    browser()
    #Pad the input catalogue with stellar contaminants down to below the magnitude limit
    dens<-density(inputmags,bw=0.5,kernel='rect',na.rm=TRUE)
    mag.mode<-dens$x[which.max(dens$y)]
    #Determine Number Counts using Driver et al r-band relation
    ngal<-length(which(inputmags<=(mag.mode+0.25) & inputmags>=(mag.mode-0.25)))/(diff(range(ra_g))*diff(range(dec_g)))
    NDriver<-10^(0.38*(mag.mode+col.corr)-4.78)
    norm<-ngal/NDriver #Normalise to represent correct area
    #Magnitude Bin Values
    x<-mag.mode+seq(-0.25,3.25,by=0.5)
    #Number of gals in each Bin
    ndraws<-norm*10^(0.38*(x+col.corr)-4.78)
    if (confuse) {
      #If we want confusion, increase this number so that
      #source density ~= psf size
      #Source Density / deg^2
      sourceDens<-(sum(ndraws)+length(inputmags))/(diff(range(ra_g))*diff(range(dec_g)))
      #Get source Density per PSF
      sDensPSF<-sourceDens*((gauss_fwhm_as/3600)^2)
      #If source density is less than probability of finding gal withing FWHM:
      if (sDensPSF<1.5) {
        ndraws<-ndraws*(1.5/sDensPSF)
      }
    }
    #Do Magnitude Draws
    y<-foreach(draws=ndraws,lo=x-0.25,hi=x+0.25,.combine='c',.export='inputmags',.inorder=FALSE)%dopar%{
      n=draws-length(which(inputmags>lo & inputmags<hi))
      if (n<=0) {
        NULL
      } else {
        runif(draws-length(which(inputmags>lo & inputmags<hi)),min=lo,max=hi)
      }
    }
    #y<-NULL
    #for( i in 1:length(x)) {
    #  draws=ndraws[i]
    #  lo=x[i]-0.25
    #  hi=x[i]+0.25
    #  y<-c(y,runif(draws-length(which(inputmags>lo & inputmags<hi)),min=lo,max=hi))
    #}
    mag.padgals<-y
    ra.padgals<-runif(length(mag.padgals),min=min(ra_g),max=max(ra_g))
    dec.padgals<-runif(length(mag.padgals),min=min(dec_g),max=max(dec_g))
    a_g.padgals<-a_g[runif(length(mag.padgals),min=1,max=length(a_g))]
    axrat.padgals<-runif(length(mag.padgals),min=min(axrat,na.rm=TRUE),max=max(axrat,na.rm=TRUE))
    b_g.padgals<-a_g.padgals*axrat.padgals
    Reff.padgals<-(a_g.padgals)*(1/(2.5*1.19))/(asperpix)
    theta.padgals<-runif(length(mag.padgals),min=min(theta_g),max=max(theta_g))
    #Add padding galaxies to the catalogue
    id_g<-c(id_g,9999000+1:length(mag.padgals))
    ra_g<-c(ra_g,ra.padgals)
    dec_g<-c(dec_g,dec.padgals)
    a_g<-c(a_g,a_g.padgals)
    b_g<-c(b_g,b_g.padgals)
    gamapos<-ad2xy(ra_g,dec_g,astr_struc)
    x_g<-gamapos[,1]
    y_g<-gamapos[,2]
    axrat<-c(axrat,axrat.padgals)
    theta_g<-c(theta_g,theta.padgals)
    Reff_pix<-c(Reff_pix,Reff.padgals)
    inputmags<-c(inputmags,mag.padgals)
    fluxweight.new<-10^c((8.9-inputmags)/2.5)
    if (any((zapsmall(fluxweight)-zapsmall(fluxweight.new[1:length(fluxweight)]))!=0,na.rm=TRUE)) {
      cat("\nFLUXWEIGHTS ARE NOT EQUAL\n")
    } else {
      fluxweight<-fluxweight.new
    }
  }

  #If needed, correct angluar coordinates
  #Details
  #Exponential Generation function uses N90E0 angular inputs
  ### This is OPPOSITE the aperture function used in LAMBDAR
  #If the input angles are not N90E0, correct for this offset
  if (!angoffset) { theta_off<-90-theta_g } else {theta_off<-theta_g}
  #Correct for any reversal of fits image
  if (astr_struc$CD[1,1]*astr_struc$CD[2,2]<0){ theta_off<-theta_off*-1 }
  #

  #ind<-order(Reff_pix)
  #Reff_pix<-Reff_pix[ind]
  #theta_off<-theta_off[ind]
  #inputmags<-inputmags[ind]
  #axrat<-axrat[ind]
  #x_g<-x_g[ind]
  #y_g<-y_g[ind]

  #Create Aperture Masks
  message("Creating Aperture Masks")
  es_mask<-foreach(mag=inputmags, theta=theta_off,axr=axrat,Reff=Reff_pix,xdelt=(x_g%%1),ydelt=(y_g%%1), .export=c("psffwhm.pix"), .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
  #es_mask<-NULL
  #for (i in 1:length(inputmags)) {
  #    mag=inputmags[i]
  #    theta=theta_off[i]
  #    axr=axrat[i]
  #    Reff=Reff_pix[i]
  #    xdelt=(x_g[i]%%1)
  #    ydelt=(y_g[i]%%1)

      #Generate Profile with Shot noise & Convolution
      if (Reff!=0) {
        #Galaxy
        #Calculate Required number of Photons
        #Formula for number of photons given magnitude, exposure time (s),
        #telescope area (m^2), and filter Effective Width and Wavelength
        photCount<-function(abmag,exp,area,Weff,lamEff) { return=10^((8.9-abmag)/2.5+7)*(1/1.51)*exp*area*(Weff/lamEff) }
        #
        #Calculate N Photons given magnitude and SDSS r-band stats
        bn=1.678
        #Get N photons
        N=photCount(mag,ObsParm$exp,ObsParm$area,ObsParm$Weff,ObsParm$lamEff)

        #Generate Image Data

        #Get N random photons from desired exponential profile
        tempr<-rexp(N, rate=bn/Reff)
        tempang<-runif(N, min=0, max=2*pi)

        #Convert from photon arrivals from Sperical to Cartesian Coords
        tempxy<-sph2car(long=tempang, lat=0, radius=tempr, deg=FALSE)

        #Stretch y-axis by ellipticity of profile
        tempxy[,2]<-(tempxy[,2])*axr

        #Rotate to desired theta
        tempxy<-rotdata2d(tempxy[,1],tempxy[,2],theta)

        #Convolve exponential profile with the PSF using bkde2D
        range.x <- list(0, 0)
        for (id in (1L:2L)) range.x[[id]] <- c(min(tempxy[, id]) - 9 * psffwhm.pix, max(tempxy[, id]) + 9 * psffwhm.pix)
        im<-bkde2D(tempxy, bandwidth=rep(psffwhm.pix/2,2), gridsize = ceiling(c(diff(range.x[[1]]),diff(range.x[[2]]))), range.x=range.x,truncate = TRUE)
        xval<-floor(im$x1)
        yval<-floor(im$x2)
        im<-zapsmall(im$fhat)
        #Make stamp square
        centre<-as.numeric(which(im==max(im), arr.ind=TRUE))
        maxim<-max(abs(c(xval,yval)))
        slen<-maxim*2+1
        delta<-floor(slen/2)*c(-1,+1)
        lims<-rbind(centre[1]+delta,centre[2]+delta)
        aplims<-rbind(c(1,slen),c(1,slen))
        #Check for -ve indexes or indexes above PSF width {{{
        if (lims[1,1]<1) {
          aplims[1,1]<-aplims[1,1]+(1-lims[1,1])
          lims[1,1]<-1
        }
        if (lims[2,1]<1) {
          aplims[2,1]<-aplims[2,1]+(1-lims[2,1])
          lims[2,1]<-1
        }
        if (lims[1,2]>length(im[,1])) {
          aplims[1,2]<-(slen-(lims[1,2]-length(im[,1])))
          lims[1,2]<-length(im[,1])
        }
        if (lims[2,2]>length(im[1,])) {
          aplims[2,2]<-(slen-(lims[2,2]-length(im[1,])))
          lims[2,2]<-length(im[1,])
        }
        #}}}
        ap<-matrix(0, nrow=slen, ncol=slen)
        ap[aplims[1]:aplims[3],aplims[2]:aplims[4]]<-im[lims[1]:lims[3],lims[2]:lims[4]]
        im<-ap
        xval<-1:length(im[,1])-ceiling(length(im[,1])/2)
        yval<-xval
      } else {
       #{{{
        #Point Source; place down PSF
        centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
        slen=length(psf[,1])
        delta<-floor(slen/2)*c(-1,+1)
        lims<-rbind(centre[1]+delta,centre[2]+delta)
        aplims<-rbind(c(1,slen),c(1,slen))
        #Check for -ve indexes or indexes above PSF width {{{
        if (lims[1,1]<1) {
          aplims[1,1]<-aplims[1,1]+(1-lims[1,1])
          lims[1,1]<-1
        }
        if (lims[2,1]<1) {
          aplims[2,1]<-aplims[2,1]+(1-lims[2,1])
          lims[2,1]<-1
        }
        if (lims[1,2]>length(psf[,1])) {
          aplims[1,2]<-(slen-(lims[1,2]-length(psf[1,])))
          lims[1,2]<-length(psf[1,])
        }
        if (lims[2,2]>length(psf[,1])) {
          aplims[2,2]<-(slen-(lims[2,2]-length(psf[,1])))
          lims[2,2]<-length(psf[,1])
        }
        #}}}
        im<-matrix(0, nrow=slen, ncol=slen)
        im[aplims[1]:aplims[3],aplims[2]:aplims[4]]<-psf[lims[1]:lims[3],lims[2]:lims[4]]
        xval<-1:length(im[,1])-ceiling(length(im[,1])/2)
        yval<-xval
        # }}}
        #range.x <- list(0, 0)
        #tempxy<-cbind(x=0,y=0,z=0)
        #for (id in (1L:2L)) range.x[[id]] <- c(min(tempxy[, id]) - 9 * psffwhm.pix, max(tempxy[, id]) + 9 * psffwhm.pix)
        #im<-bkde2D(tempxy, bandwidth=rep(psffwhm.pix/2,2), gridsize = ceiling(c(diff(range.x[[1]]),diff(range.x[[2]]))), range.x=range.x,truncate = TRUE)
        #xval<-floor(im$x1)
        #yval<-floor(im$x2)
        #im<-zapsmall(im$fhat)
      }
      #Scale to desired Input Magnitude
      im<-im/sum(im)*10^((8.9-mag)/2.5)

      #
      #return=im
      #If needed, Interpolate onto correct grid
      if ((length(im)>1)) { # &((xdelt!=0.5)|(ydelt!=0.5))) {
        len<-dim(im)
        #Make grid at old pixel centres
        old_obj<-list(x=seq(1,len[1]), y=seq(1,len[2]),z=im)
        #Make expanded grid of new pixel centres
        expanded<-expand.grid(seq(1,len[1]),seq(1,len[2]))
        xnew<-expanded[,1]-xdelt
        ynew<-expanded[,2]-ydelt
        #
        #Interpolate
        im<-matrix(interp2D(xnew, ynew, old_obj), ncol=len[2])
        #
        #
      }
      return=im
      #es_mask<-c(es_mask, list(im))
  }
  message("Aperture Creation Complete")
  #

  #Check for Bad objects after generation
  es_mask<-foreach(smask=es_mask, .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
    #Check masking to determine if Aperture is acceptable
    if (any(is.na(smask))) {
      warning("Mask contains NA values. Setting those pixels to 0")
      smask[which(is.na(smask))]<-0
    }
    smask
  }
  stamplen<-foreach(smask=es_mask, .inorder=TRUE, .options.mpi=mpiopts, .combine='c') %dopar% { length(smask[,1]) }

  #Check that apertures do not cross image mask boundary
  if (((length(image.env$imm)!=1)&(length(which(image.env$imm!=1))!=0))) {
    #Check Mask stamps for Aperture Acceptance
    message('Combining Aps with Mask Stamps')
    esa_mask<-foreach(slen=stamplen, smask=es_mask,mmask=imm_mask, .export="useMaskLim", .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
      #Check masking to determine if Aperture is acceptable
      if (any(is.na(smask))) {
        warning("Mask contains NA values. Setting those pixels to 0")
        smask[which(is.na(smask))]<-0
      }
      if (mean(mmask)<useMaskLim) {
        #Too much. Skip
        array(0, dim=c(slen,slen))
        #
      } else {
        #Not too much. Keep
        smask
        #
      }
      #
    }
    message('Mask Combine Finished.')
    #
  } else {
    #No image mask, all can be kept
    esa_mask<-es_mask
    #
  }#

  #Determine Stamp Limits in Image Space
  factor<-rep(0,length(esa_mask))
  for (i in 1:length(esa_mask)) { factor[i]<-((length(esa_mask[[i]][,1])-1)/2) }
  image_lims<-cbind(floor(x_g)-factor, floor(x_g)+factor, floor(y_g)-factor, floor(y_g)+factor)
  ind<-which(image_lims[,1]<1)
  for (i in ind) {
    esa_mask[[i]]<-esa_mask[[i]][(image_lims[i,1]-1):0,]
    image_lims[i,1]<-1
  }
  ind<-which(image_lims[,3]<1)
  for (i in ind) {
    esa_mask[[i]]<-esa_mask[[i]][,(image_lims[i,3]-1):0]
    image_lims[i,3]<-1
  }
  ind<-which(image_lims[,2]>length(image.env$im[,1]))
  for (i in ind) {
    over<-image_lims[i,2]-length(image.env$im[,1])
    lim<-length(esa_mask[[i]][,1])-over
    esa_mask[[i]]<-esa_mask[[i]][1:lim,]
    image_lims[i,2]<-length(image.env$im[,1])
  }
  ind<-which(image_lims[,4]>length(image.env$im[1,]))
  for (i in ind) {
    over<-image_lims[i,4]-length(image.env$im[1,])
    lim<-length(esa_mask[[i]][1,])-over
    esa_mask[[i]]<-esa_mask[[i]][,1:lim]
    image_lims[i,4]<-length(image.env$im[1,])
  }
  #

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(esa_mask)))/length(esa_mask)*100,digits=2),"% of the esa_mask matrix are NA"))
  } #

  #Parse Parameter Space
  if (!is.null(env)) { detatch(env) }
  assign("id_g",id_g,envir=outenv)
  assign("ra_g",ra_g,envir=outenv)
  assign("dec_g",dec_g,envir=outenv)
  assign("a_g",a_g,envir=outenv)
  assign("b_g",b_g,envir=outenv)
  assign("x_g",x_g,envir=outenv)
  assign("y_g",y_g,envir=outenv)
  assign("theta_g",theta_g,envir=outenv)
  assign("Reff_pix",Reff_pix,envir=outenv)
  assign("axrat",axrat,envir=outenv)
  assign("inputmags",inputmags,envir=outenv)
  assign("fluxweight",fluxweight,envir=outenv)
  assign("image_lims",image_lims,envir=outenv)
  #

  message('===========END============Make_ESA_MASK=============END=================\n')

  #Return array of apertures
  return=esa_mask
  #
}

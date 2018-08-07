make.exponential.apertures <-
function(outenv=parent.env(environment()), env=NULL,ObsParm,padGals,col.corr=0,confuse=FALSE,subs=NULL){
#Procedure creates exponential galaxy profiles, using input parameters
#from the catalogue, and places them (in order) onto stamps
#Procedure is parallelised to allow scaleability

  message("Creating Simulated Galaxies from Input Catalogue")
  message("Assumptions are:\nCatalogue Galaxy Apertures are 2.5*Kron Apertures\nCatalogue Fluxweights are Jy fluxes for use in Simulation")
  if (!quiet) { cat('Make_ESA_Mask   ') }
  message('--------------------------Make_ESA_Mask---------------------------------')

  if (!is.null(subs)) { 
    subs<-1:length(cat.id)
  } 
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
  Reff_pix<-(cat.a)*(1/(2.5*kronrad(1)))/(arcsec.per.pix)
  #

  #Aperture Axis Ratio in [0,1]
  axrat<-cat.b/cat.a
  axrat[which(!is.finite(axrat))]<-1
  #

  #Get input magnitudes
  if (is.na(mag.zp)) {
    stop("No Magnitude Zero Point supplied or read. magnitudes are required for Sim")
  }
  inputmags<- -2.5*log10(flux.weight)+mag.zp
  #

  #Remove Galaxies with non-finite input magnitudes
  subs<-subs[which(is.finite(inputmags[subs]))]
  #inputmags<-inputmags[index]
  #cat.a<-cat.a[index]
  #cat.b<-cat.b[index]
  #cat.ra<-cat.ra[index]
  #cat.dec<-cat.dec[index]
  #cat.x<-cat.x[index]
  #cat.y<-cat.y[index]
  #cat.theta<-cat.theta[index]
  #Reff_pix<-Reff_pix[index]
  #flux.weight<-flux.weight[index]
  #axrat<-axrat[index]
  #cat.id<-cat.id[index]
  #contams<-contams[index]
  #

  if (padGals) {
    #Pad the input catalogue with stellar contaminants down to below the magnitude limit
    dens<-density(inputmags,bw=0.1,na.rm=TRUE,from=0,to=50)
    mag.mode<-dens$x[which.max(dens$y)]
    #Determine Number Counts using Driver et al r-band relation
    x<-seq(quantile(inputmags,c(0.01)),mag.mode+2.25,by=0.5)
    NDriver<-10^(-0.01476*x^2+1.025*x-11.76)
    #Normalise for sky area
    skylims<-dim(image.env$im)*arcsec.per.pix/3600
    skyarea<-skylims[1]*skylims[2]
    #Magnitude Bin Values
    NDriver<-skyarea*NDriver
    #Number of gals in each Bin
    ndraws<-rpois(length(NDriver),NDriver)
    cat("Source Density: ",(sum(ndraws))/skyarea,' Objects/sqDegree\n')
    if (confuse) {
      #If we want confusion, increase this number so that
      #source density ~= psf size
      #Source Density / deg^2
      sourceDens<-(sum(ndraws))/skyarea
      #Get source Density per PSF
      sDensPSF<-sourceDens*((gauss.fwhm.arcsec/3600)^2)
      #If source density is less than probability of finding gal withing FWHM:
      if (sDensPSF<1.5) {
        add=0
        while (sDensPSF<1.5) {
          add=add+0.5
          #Get another bin of draws
          x<-seq(quantile(inputmags,c(0.01)),mag.mode+2.25+add,by=0.5)
          #Get new draws
          NDriver<-10^(-0.01476*x^2+1.025*x-11.76)
          NDriver<-skyarea*NDriver
          ndraws<-rpois(length(NDriver),NDriver)
          #Source Density / deg^2
          sourceDens<-(sum(ndraws))/skyarea
          #Get source Density per PSF
          sDensPSF<-sourceDens*((gauss.fwhm.arcsec/3600)^2)
        }
        #ndraws<-ndraws*(1.5/sDensPSF)
      }
      cat("Confused Source Density: ",sum(ndraws)/skyarea,' Objects/sqDegree\n')
    }
    #Do Magnitude Draws
    y<-foreach(draws=ndraws,lo=x-0.25,hi=x+0.25,.combine='c',.export='inputmags',.inorder=FALSE)%dopar%{
      n=draws-length(which(inputmags>=lo & inputmags<hi))
      if (n<=0) {
        NULL
      } else {
        runif(n,min=lo,max=hi)
      }
    }
    y<-y+rnorm(length(y),0,0.2)
    photCount<-function(abmag,exp,area,Weff,lamEff) { return=10^((8.9-abmag)/2.5+7)*(1/1.51)*exp*area*(Weff/lamEff) }
    N=photCount(x,ObsParm$exp,ObsParm$area,ObsParm$Weff,ObsParm$lamEff)
    cat("Photon Counts in each x-bin:\nBin Centre: ",x,"\nN Phot: ",N,"\n")
    mag.padgals<-y
    ra.padgals<-runif(length(mag.padgals),min=min(cat.ra),max=max(cat.ra))
    dec.padgals<-runif(length(mag.padgals),min=min(cat.dec),max=max(cat.dec))
    cat.a.padgals<-cat.a[runif(length(mag.padgals),min=1,max=length(cat.a))]
    axrat.padgals<-runif(length(mag.padgals),min=min(axrat,na.rm=TRUE),max=1)
    cat.b.padgals<-cat.a.padgals*axrat.padgals
    Reff.padgals<-(cat.a.padgals)*(1/(2.5*kronrad(1)))/(arcsec.per.pix)
    theta.padgals<-runif(length(mag.padgals),min=min(cat.theta),max=max(cat.theta))
    #Add padding galaxies to the catalogue
    cat.id<-c(cat.id,9999000+1:length(mag.padgals))
    cat.ra<-c(cat.ra,ra.padgals)
    cat.dec<-c(cat.dec,dec.padgals)
    cat.a<-c(cat.a,cat.a.padgals)
    cat.b<-c(cat.b,cat.b.padgals)
    gama.pos<-ad.to.xy(cat.ra,cat.dec,astr.struc)
    cat.x<-gama.pos[,1]
    cat.y<-gama.pos[,2]
    axrat<-c(axrat,axrat.padgals)
    cat.theta<-c(cat.theta,theta.padgals)
    Reff_pix<-c(Reff_pix,Reff.padgals)
    inputmags<-c(inputmags,mag.padgals)
    flux.weight<-10^c((8.9-inputmags)/2.5)
    contams<-c(contams,rep(1,length(mag.padgals)))
  }

  #If needed, correct angluar coordinates
  #Details
  #Exponential Generation function uses N90E0 angular inputs
  ### This is OPPOSITE the aperture function used in LAMBDAR
  #If the input angles are not N90E0, correct for this offset
  if (!ang.offset) { theta.offset<-90-cat.theta } else {theta.offset<-cat.theta}
  #Correct for any reversal of fits image
  if (astr.struc$CD[1,1]*astr.struc$CD[2,2]<0){ theta.offset<-theta.offset*-1 }

  print(length(inputmags))

  ##Reorder for faster parallelisation
  #nchild<-num.cores
  #ind<-order(inputmags)
  #matdim<-c(ceiling(length(inputmags)/nchild),nchild)
  #ind<-c(ind,rep(NA,matdim[1]*matdim[2]-length(ind)))
  #ind<-matrix(ind,ncol=matdim[2],nrow=matdim[1],byrow=TRUE)
  #ind<-as.numeric(ind)
  #ind<-ind[which(!is.na(ind))]

  #Reff_pix<-Reff_pix[ind]
  #theta.offset<-theta.offset[ind]
  #inputmags<-inputmags[ind]
  #axrat<-axrat[ind]
  #contams<-contams[ind]
  #cat.ra<-cat.ra[ind]
  #cat.dec<-cat.dec[ind]
  #cat.x<-cat.x[ind]
  #cat.y<-cat.y[ind]
  #cat.id<-cat.id[ind]
  #cat.a<-cat.a[ind]
  #cat.b<-cat.b[ind]
  #cat.theta<-cat.theta[ind]
  #flux.weight<-flux.weight[ind]

  #Create Aperture Masks
  message("Creating Aperture Masks")
  psfsigma.pix<-psffwhm.pix/(2*sqrt(2*log(2)))
  if (confidence==1) {
    warning("Cannot have PSF Confidence = Unity when analytically deriving PSF; Using Maximum double value (1-1E-16)")
    confidence=1-1E-16
  }
  nsig<-sqrt(qchisq(confidence,df=2)) #Determine nsigma for desired confidence using chisq distribution
  psf.clip<-ceiling(nsig*psfsigma.pix)
  psf.clip<-psf.clip*2+1 # convert psf.clip from radius to diameter (and make sure it's odd)
  es_mask<-rep(list(NULL),length(inputmags))
  pb<-txtProgressBar(min=1,max=length(inputmags),style=3)
  for (i in subs) {
      setTxtProgressBar(pb,i)
      mag=inputmags[i]
      theta=theta.offset[i]
      axr=axrat[i]
      Reff=Reff_pix[i]
      xdelt=(cat.x[i]%%1)
      ydelt=(cat.y[i]%%1)

      #Calculate Required number of Photons
      #Formula for number of photons given magnitude, exposure time (s),
      #telescope area (m^2), and filter Effective Width and Wavelength
      photCount<-function(abmag,exp,area,Weff,lamEff) { return=10^((8.9-abmag)/2.5+7)*(1/1.51)*exp*area*(Weff/lamEff) }
      #
      #Calculate N Photons given magnitude and SDSS r-band stats
      bn=1.678
      #Get N photons
      N=floor(photCount(mag,ObsParm$exp,ObsParm$area,ObsParm$Weff,ObsParm$lamEff))
      if (N>1E4) {
      #Use Analytic {{{
        stamplen<-(floor((ceiling(Reff*9)+ceiling(psf.clip))/2)*2+5)
        xy<-expand.grid(1:(stamplen*2),1:(stamplen*2))
        tempxy<-xy
        xy[,1]<-xy[,1]-stamplen-0.5
        xy[,2]<-xy[,2]-stamplen-0.5
        tempxy[,1]<-tempxy[,1]-stamplen-0.5-xdelt
        tempxy[,2]<-tempxy[,2]-stamplen-0.5-ydelt
        if (Reff!=0) {
          #Galaxy
          tempxy<-rotate.data.2d(tempxy[,1],tempxy[,2],90-theta)
          tempxy[,1]<-tempxy[,1]/axr
          I<-exp(-bn*(sqrt((xy[,1])^2+(xy[,2])^2)/Reff-1))
          I[which(I<max(I)*1E-5)]<-0
          im<-zapsmall(matrix(interp.2d(tempxy[,1],tempxy[,2], list(x=1:(stamplen*2)-stamplen-0.0,y=1:(stamplen*2)-stamplen-0.0,z=matrix(zapsmall(I),stamplen*2,stamplen*2)))[,3], ncol=stamplen*2,nrow=stamplen*2,byrow=FALSE))
          im<-im[(ceiling(stamplen*0.5)):(floor(stamplen*1.5)),(ceiling(stamplen*0.5)):(floor(stamplen*1.5))]
          #Convert to Jy
          im<-im/sum(im)*10^((8.9-mag)/2.5)
        }
        if (psf.filt || Reff==0) {
        #Point Source
          ps<-matrix(dnorm(sqrt((xy[,1])^2+(xy[,2])^2),mean=0,sd=psfsigma.pix),stamplen*2,stamplen*2)
          ps<-ps[(ceiling(stamplen*0.5)):(floor(stamplen*1.5)),(ceiling(stamplen*0.5)):(floor(stamplen*1.5))]
          if (Reff==0) {
            #Convert to Jy
            im<-ps/sum(ps)*10^((8.9-mag)/2.5)
          } else {
            #Convert to Jy
            im<-convolve.psf(ps,im)
            im<-im/sum(im)*10^((8.9-mag)/2.5)
          }
        }
        es_mask[[i]]<-im
      #}}}
      } else if (N > 0) {
        #Generate Profile with Shot noise & Convolution
        if (Reff!=0) {
          #Galaxy
          #Get N random photons from desired exponential profile
          tempr<-rexp(N, rate=bn/Reff)
          tempang<-runif(N, min=0, max=2*pi)

          #Convert from photon arrivals from Sperical to Cartesian Coords
          tempxy<-sph.to.car(long=tempang, lat=0, radius=tempr, deg=FALSE)

          #Stretch y-axis by ellipticity of profile
          tempxy[,2]<-(tempxy[,2])*axr

          #Rotate to desired theta
          tempxy<-rotate.data.2d(tempxy[,1],tempxy[,2],theta)

          #Filter through Gaussian PSF
          if (psf.filt) {
            tempxy[1:length(tempxy)]<-tempxy[1:length(tempxy)]+rnorm(length(tempxy),mean=0,sd=psfsigma.pix)
          }
        } else {
          #Point Source; Use Gaussian with PSF FWHM
          tempxy<-matrix(rnorm(2*N,mean=0,sd=psfsigma.pix),ncol=2)
        }
        #Centroid correctly
        stamplen<-(floor((ceiling(Reff*9)+ceiling(psf.clip))/2)*2+5)
        tempx<-(tempxy[,1]+xdelt+stamplen/2)+0.5
        tempy<-(tempxy[,2]+ydelt+stamplen/2)+0.5
        k<-kde2d(tempx,tempy,lims=c(1,stamplen,1,stamplen),n=stamplen)
        im<-matrix(k$z,stamplen,stamplen)

        #Convert from ADU to Jy
        im<-im/sum(im)*10^((8.9-mag)/2.5)
        es_mask[[i]]<-im
      } 
  }
  close(pb)
  message("Aperture Creation Complete")
  #

  #Check for Bad objects after generation
  es_mask<-foreach(smask=es_mask, .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
    #Check masking to determine if Aperture is acceptable
    if (any(is.na(smask))) {
      warning("Mask contains NA values. Setting those pixels to 0")
      smask[which(is.na(smask))]<-0
    }
    smask
  }
  stamplen<-foreach(smask=es_mask, .inorder=TRUE, .options.mpi=mpi.opts, .combine='c') %dopar% { length(smask[,1]) }

  #Check that apertures do not cross image mask boundary
  if (((length(image.env$imm)!=1)&(length(which(image.env$imm!=1))!=0))) {
    #Check Mask stamps for Aperture Acceptance
    message('Combining Aps with Mask Stamps')
    esa_mask<-foreach(slen=stamplen, smask=es_mask,mmask=mask.stamp, .export="use.mask.lim", .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
      #Check masking to determine if Aperture is acceptable
      if (any(is.na(smask))) {
        warning("Mask contains NA values. Setting those pixels to 0")
        smask[which(is.na(smask))]<-0
      }
      if (mean(mmask)<use.mask.lim) {
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
  ap.lims.data.map<-cbind(floor(cat.x)-factor, floor(cat.x)+factor, floor(cat.y)-factor, floor(cat.y)+factor)
  ind<-which(ap.lims.data.map[,1]<1)
  for (i in ind) {
    esa_mask[[i]]<-esa_mask[[i]][(ap.lims.data.map[i,1]-1):0,]
    ap.lims.data.map[i,1]<-1
  }
  ind<-which(ap.lims.data.map[,3]<1)
  for (i in ind) {
    esa_mask[[i]]<-esa_mask[[i]][,(ap.lims.data.map[i,3]-1):0]
    ap.lims.data.map[i,3]<-1
  }
  ind<-which(ap.lims.data.map[,2]>length(image.env$im[,1]))
  for (i in ind) {
    over<-ap.lims.data.map[i,2]-length(image.env$im[,1])
    lim<-length(esa_mask[[i]][,1])-over
    esa_mask[[i]]<-esa_mask[[i]][1:lim,]
    ap.lims.data.map[i,2]<-length(image.env$im[,1])
  }
  ind<-which(ap.lims.data.map[,4]>length(image.env$im[1,]))
  for (i in ind) {
    over<-ap.lims.data.map[i,4]-length(image.env$im[1,])
    lim<-length(esa_mask[[i]][1,])-over
    esa_mask[[i]]<-esa_mask[[i]][,1:lim]
    ap.lims.data.map[i,4]<-length(image.env$im[1,])
  }
  #

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(esa_mask)))/length(esa_mask)*100,digits=2),"% of the esa_mask matrix are NA"))
  } #

  #Parse Parameter Space
  if (!is.null(env)) { detach(env) }
  if (!is.null(subs)) { 
  } 
  assign("cat.id",cat.id,envir=outenv)
  assign("cat.ra",cat.ra,envir=outenv)
  assign("cat.dec",cat.dec,envir=outenv)
  assign("cat.a",cat.a,envir=outenv)
  assign("cat.b",cat.b,envir=outenv)
  assign("cat.x",cat.x,envir=outenv)
  assign("cat.y",cat.y,envir=outenv)
  assign("contams",contams,envir=outenv)
  assign("cat.theta",cat.theta,envir=outenv)
  assign("Reff_pix",Reff_pix,envir=outenv)
  assign("axrat",axrat,envir=outenv)
  assign("inputmags",inputmags,envir=outenv)
  assign("flux.weight",flux.weight,envir=outenv)
  assign("ap.lims.data.map",ap.lims.data.map,envir=outenv)

  message('===========END============Make_ESA_MASK=============END=================\n')

  #Return array of apertures
  return=esa_mask
  #
}

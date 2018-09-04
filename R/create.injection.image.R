#
#Create Simulated Image
#
create.injection.image<-
function(n.inject=100, ObsParm, col.corr=0, outenv=parent.env(environment()), env=NULL){

  #Setup Environments {{{
  environment(make.exponential.apertures)<-environment()
  environment(make.aperture.map)<-environment()
  environment(read.psf)<-environment()
  for (nam in ls.deb("package:LAMBDAR",simple=TRUE)) { 
    if (nam%in%ls(envir=environment())) { 
      debug(get(nam,envir=environment()))
    }
  }
  #}}}

  #Get PSF details {{{
  timer<-proc.time()
  if (!quiet) { cat('Getting PSF details') }
  message('Getting PSF details')
  #Details {{{
  #Calculate PSF - if one has not be supplied,
  #then a gaussian PSF will be created. Stampsize of PSF
  #should be = maximum stampsize: max aperture * stamp mult }}}
  if (gauss.fwhm.arcsec==0) {
    gauss.fwhm.arcsec=get.fwhm(read.fits(file.path(path.root,path.work,psf.map))$dat[[1]])*arcsec.per.pix
  }
  #Get PSF {{{
  psf<-read.psf(outenv=environment(),"NONE",arcsec.per.pix,max(cat.a,na.rm=TRUE),1-(1E-14),gauss.fwhm.arcsec=gauss.fwhm.arcsec)
  if (sim.gauss.arcsec>0) {
    sim.gauss.psf<-read.psf(outenv=environment(),"NONE",arcsec.per.pix,max(cat.a,na.rm=TRUE),1-(1E-14),gauss.fwhm.arcsec=sim.gauss.arcsec)
  } else {
    sim.gauss.psf<-0
  }
  #}}}
  #Get radius of FWHM using FWHM confidence value = erf(2*sqrt(2*log(2))/sqrt(2)) {{{
  psffwhm.pix<-get.fwhm(psf)
  #}}}
  #}}}
  #Notify {{{
  if (showtime) { cat(paste(' - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )\n'))
    message(paste('Getting PSF details - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat(' - Done\n') }
  #}}}

  #Get Input Image Noise Characteristics {{{
  timer<-proc.time()
  if (verbose) { cat("Determine Image Noise Characteristics ") }
  if (length(image.env$im) > 1E6) { index<-runif(1E6,min=1,max=length(image.env$im)) } else { index<-1:length(image.env$im) }
  xdat<-as.numeric(image.env$im[index])
  ind<-which(xdat < mean(xdat)*1.1)
  if (length(ind)>2) {
    x.dens=density(xdat[ind])
    x.mode=x.dens$x[which.max(x.dens$y)]
    # Select all data with value < mode - this will be exclusively noise
    xdat=xdat-x.mode
    noise=xdat[which(xdat <= 0.0)]
    # Determine Gaussian Characteristics
    stdev=sd(noise)
    message("Noise Profile: Mean = ", x.mode, " ; Stdv = ",abs(stdev)," ;\n")
  } else {
    stop("Badness in Image Noise characteristic determination")
  }#}}}
  #Notify {{{
  if (showtime) { cat(paste(' - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )\n'))
  } else if (!quiet) { cat(' - Done\n') }
  #}}}

  #Create the injection sources {{{
  contams<-c(rep(TRUE,length(cat.x)),rep(FALSE,n.inject))
  cat.id<-c(cat.id,paste0("Injected_",1:n.inject))
  samp.x<-sample(1:length(cat.x),n.inject)
  samp.y<-sample(1:length(cat.y),n.inject)
  cat.x<-c(cat.x,cat.x[samp.x])
  cat.y<-c(cat.y,cat.y[samp.y])
  cat.ra<-c(cat.ra,cat.ra[samp.x])
  cat.dec<-c(cat.dec,cat.dec[samp.y])
  cat.theta<-c(cat.theta,runif(n.inject,min=-90,max=90))
  cat.a<-c(cat.a,runif(n.inject,min=0,max=quantile(cat.a,0.9,na.rm=T)))
  cat.b<-c(cat.b,runif(n.inject,min=0.2,1)*cat.a[which(!contams)])
  cat("stdev:", stdev,'\n')
  cat("med flux:", median(runif(n.inject,min=5,max=100)*stdev*10*pi*(cat.a*cat.b)[which(!contams)]),'\n')
  sim.fluxes<-runif(n.inject,min=5,max=100)*stdev*10*pi*(cat.a*cat.b)[which(!contams)]
  flux.true<-c(flux.weight,sim.fluxes)
  flux.weight<-c(rep(min(flux.weight,na.rm=T)/1000,length(flux.weight)),sim.fluxes)
  flux.weight<-flux.true
  contams<-rep(FALSE,length(contams))
  #}}}

  #Create Simulated Profiles & Image {{{
  if (!quiet) { cat(paste('Creating Injection Image  ')) }
  timer=system.time(esa<-make.exponential.apertures(outenv=environment(),ObsParm=ObsParm,padGals=FALSE,col.corr=col.corr,confuse=FALSE))
  simFlux<-foreach(esam=esa, .inorder=TRUE, .options.mpi=mpi.opts, .combine='c', .noexport=ls(envir=environment())) %dopar% { sum(esam) }
  simFlux[which(contams)]<-flux.true[which(contams)]
  npix<-foreach(esam=esa, .inorder=TRUE, .options.mpi=mpi.opts, .combine='c', .noexport=ls(envir=environment())) %dopar% { length(esam) }
  cat.out<-data.frame(CATAID=cat.id,RA=cat.ra,DEC=cat.dec,X.im=cat.x,Y.im=cat.y,REff=Reff_pix*arcsec.per.pix,SEMIMAJ.arcsec=cat.a,
           SEMIMIN.arcsec=cat.b,THETA=cat.theta,inputFlux=flux.weight,StampFlux=simFlux,nPix=npix,CONTAM=ifelse(contams,1,0))
  colnames(cat.out)<-c(cata.lab,ra.lab,dec.lab,"X.im","Y.im","R.Eff",semimaj.lab,semimin.lab,theta.lab,"Input.Flux.Jy","Sim.Flux.Jy","n.pix",contam.lab)
  write.csv(file=file.path(path.root,path.work,path.out,"SimFlux.csv"), cat.out, row.names=FALSE,quote=FALSE)
  timer1=system.time(image.env$ea<-make.aperture.map(outenv=environment(), esa, dim(image.env$im),subs=which(!contams)))
  if (showtime) { cat("   - Done (",round(timer[3]+timer1[3],digits=2),"sec )\n")
    message(paste('Create Sim Image - Done (',round(timer[3]+timer1[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #Convert Noise mode & stdev to Jys
  message("Converting Sim Image to Jy from ADU.")
  x.mode<-x.mode*10^((8.9-mag.zp)/2.5)
  stdev<-stdev*10^((8.9-mag.zp)/2.5)
  image.env$data.hdr[gain.label,1]<-gain*10^((8.9-mag.zp)/2.5)
  image.env$data.hdr[satur.label,1]<-saturation*10^((8.9-mag.zp)/2.5)
  message("Noise Properties in input Image (Jy): mode=",x.mode,"; sd=",stdev)
  if (!quiet) { cat("   - Done\n") }
  if (!quiet) { cat(paste('Outputting Simulate Image to',"sim_image.fits","   ")) }
  timer=system.time(write.fits.image.file(file.path(path.root,path.work,path.out,"sim_image.fits"),(image.env$ea+0*image.env$im*10^((8.9-mag.zp)/2.5)),image.env$data.hdr,nochange=TRUE) )
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Output Sim Image - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}

  return=NULL

}

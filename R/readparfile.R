readparfile <-
function(parfile=NA, starttime=NA, quiet=FALSE, env=NULL){
  #Procedure to setup Parameter Space (Read .par File)
  if (is.null(env)) { env<-environment() }

  #Check Calling Syntax
  if (is.na(parfile)) {
    stop("Parameter file not supplied. To create the default parameter file, run MeasureGamaFluxes('--makepar').")
  }
  if (is.na(starttime)) {
    warning("Start time not supplied - using current clock time")
    starttime=proc.time()[3]
  }

  #------------------------------------------------------------------------------------------------

  #Print Header
  if (!quiet) { cat(paste('------------------------------------------------------\n'))
                cat(paste('   LAMBDAR     version ',packageVersion("LAMBDAR"),': MeasureGamaFluxes\n'))
                cat(paste('------------------------------------------------------\n'))
                cat('Initialising Workspace {\n')
                cat('   Reading Parameter File   ') }
  #Test Reading of Parameter File
  error<-try(read.table(parfile, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#", row.names=1), silent=TRUE)
  if (class(error)=="try-error") {
    #Stop on Error
    #geterrmessage()
    stop("Parameter file read failed")
  }
  #Read Parameter File
  params<-read.table(parfile, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#", row.names=1)
  if (!quiet) { cat('- Done\n   Assigning Parameter Variables') }

  ## Assign Parameter File values to variables
  #Root Directory path
  ID="RootDirectory"
  pathroot<-params[ID,1]
  if (is.na(pathroot)) {
  stop("Root Directory Path not in Parameter File")
  }
  #Ensure path ends in a '/'
  if (lastnchar(pathroot,1) != '/') { pathroot<-paste(pathroot,'/') }

  #Output Directory path
  ID="OutputDirectory"
  pathout<-params[ID,1]
  if (is.na(pathout)) {
  stop("Output Directory Path not in Parameter File")
  }
  #Ensure poth ends in a '/'
  if (lastnchar(pathout,1) != '/') { pathout<-paste(pathout,'/') }
  dir.create(paste(pathroot, pathout,sep=""), showWarnings = FALSE)

  #Beam area in square arcsec
  ID="BeamArea_SqAS"
  beamarea_SOM_as<-as.numeric(params[ID,1])
  if (is.na(beamarea_SOM_as)) { 
    warning("Beamarea Parameter not present in Parameter File")
    beamarea_SOM_as<-0 
  }

  #Do we want to Convolve the apertures with a PSF
  ID="PSFConvolve"
  psffilt<-as.numeric(params[ID,1])
  if (is.na(psffilt)) {
    warning("PSF Convolve Flag not in Parameter File")
    psffilt<-0
  }
  if (psffilt==1) { nopsf<-0 } else { nopsf <- 1 }

  #PSF map filename
  ID="PSFMap"
  psfmap<-params[ID,1]
  if (is.na(psfmap)) {
    warning("PSF Map Filename not in Paramter File")
    psfmap<-"NONE"
  }

  if (psfmap=="NONE") {
    #FWHM of seeing gaussian
    ID="Gauss_FWHM_AS"
    gauss_fwhm_as<-as.numeric(params[ID,1])
    if (is.na(gauss_fwhm_as)) { gauss_fwhm_as<-0.0 }

    #If we want convolution, there is no PSF, and no gaussian FWHM provided - ERROR
    if ((psffilt)&(gauss_fwhm_as==0.0)&(psfmap=="NONE")) {
      stop("Parameter file does not provide either PSF map or Gaussian FWHM")
    }
  } else {
    #If there is a PSF Map - set gauss fwhm to zero
    gauss_fwhm_as<-0.0
  }

  #Output Ellipse Overlaid Image?
  overlaymap<-NULL
  ID="OverlayEllipse"
  overlay<-params[ID,1]
  if (is.na(overlay)) {
    warning("Overlay-Ellipses Image Flag not in Parameter File")
    overlay<-FALSE
  } else { overlay<-(overlay==1) }

  if ( overlay ) {
    #Name of the output Residual image
    ID="OverlayFile"
    overlaymap<-params[ID,1]
    if (is.na(overlaymap)) {
      warning("Overlay Ellipse Map Filename not in Parameter File")
      overlaymap<-"OverlayEllipseImage.fits"
    }
  }

  #Perform Contaminant removal?
  nocontammap<-NULL
  ID="RemoveContam"
  filtcontam<-params[ID,1]
  if (is.na(filtcontam)) {
    warning("Remove Contaminants Flag not in Parameter File")
    filtcontam<-FALSE
  } else { filtcontam<-(filtcontam==1) }

  if ( filtcontam ) {
    #Name of the output Residual image
    ID="NoContamImageFile"
    nocontammap<-params[ID,1]
    if (is.na(nocontammap)) {
      warning("Contaminant Subtracted Residual Map Filename not in Parameter File")
      nocontammap<-"NoContamResidualImage.fits"
    }
  }

  #Name of Source Catalogue
  ID="Catalogue"
  catalogue<-params[ID,1]
  if (is.na(catalogue)) {
    stop("Catalogue Path not in Parameter File")
  }

  #What is the title of the Catalogue's RA Column?
  ID="CatIDColumnLabel"
  catalab<-params[ID,1]
  if (is.na(catalab)) {
    warning("Catalogue CATID Column Label not in Parameter File; using 'CATAID'")
    ralab<-"CATAID"
  }
  #What is the title of the Catalogue's RA Column?
  ID="RAColumnLabel"
  ralab<-params[ID,1]
  if (is.na(ralab)) {
    warning("Catalogue RA Column Label not in Parameter File; using 'ALPHA_J2000'")
    ralab<-"ALPHA_J2000"
  }
  #What is the title of the Catalogue's Dec Column?
  ID="DecColumnLabel"
  declab<-params[ID,1]
  if (is.na(declab)) {
    warning("Catalogue Dec Column Label not in Parameter File; using 'DELTA_J2000'")
    declab<-"DELTA_J2000"
  }
  #What is the title of the Catalogue's Theta Column?
  ID="ThetaColumnLabel"
  thetalab<-params[ID,1]
  if (is.na(thetalab)) {
    warning("Catalogue Theta Column Label not in Parameter File; using 'THETA_J2000'")
    thetalab<-"THETA_J2000"
  }
  #What is the title of the Catalogue's SemiMaj Axis Column?
  ID="SemiMajColumnLabel"
  semimajlab<-params[ID,1]
  if (is.na(semimajlab)) {
    warning("Catalogue SemiMajor Axis Column Label not in Parameter File; using 'SEMIMAJ_AS'")
    semimajlab<-"SEMIMAJ_AS"
  }
  #What is the title of the Catalogue's SemiMaj Axis Column?
  ID="SemiMinColumnLabel"
  semiminlab<-params[ID,1]
  if (is.na(semiminlab)) {
    warning("Catalogue SemiMinor Axis Column Label not in Parameter File; using 'SEMIMIN_AS'")
    semiminlab<-"SEMIMIN_AS"
  }
    #What is the title of the Catalogue's SemiMaj Axis Column?
  ID="ContamColumnLabel"
  contamlab<-params[ID,1]
  if (is.na(contamlab)) {
    warning("Catalogue SemiMajor Axis Column Label not in Parameter File; using 'CONTAM'")
    contamlab<-"CONTAM"
  }
    #What is the title of the Catalogue's SemiMaj Axis Column?
  ID="FluxWgtColumnLabel"
  fluxweightlab<-params[ID,1]
  if (is.na(fluxweightlab)) {
    warning("Catalogue FluxWeight Column Label not in Parameter File; using 'FLUXWEIGHT'")
    fluxweightlab<-"FLUXWEIGHT"
  }
  
  #Name of Data Image
  ID="DataMap"
  datamap<-params[ID,1]
  if (is.na(datamap)) {
    stop("Data Map Path not in Parameter File")
  }

  #Name of Error Map
  ID="ErrorMap"
  errormap<-params[ID,1]
  if (is.na(errormap)) {
    warning("Error Map Path not in Parameter File")
    errormap<-"NONE"
  }

  #Name of Mask Map
  ID="MaskMap"
  maskmap<-params[ID,1]
  if (is.na(maskmap)) {
    warning("Mask Map Path not in Parameter File")
    maskmap<-"NONE"
  }
  
  #Name of Weight Map
  ID="WeightMap"
  wgtmap<-params[ID,1]
  if (is.na(wgtmap)) {
    warning("Weight Map Path not in Parameter File")
    wgtmap<-"NONE"
  }
  
  #Zero Point of Weight Map
  ID="WeightMapZP"
  wgtzp<-as.numeric(params[ID,1])
  if (is.na(wgtzp)) {
    warning("Weight Map Zero Point not in Parameter File; Using 0")
    wgtzp<-0
  }

  #Extension number of Data in FITS Header
  ID="DataExtn"
  extn<-as.numeric(params[ID,1])
  if (is.na(extn)) {
    warning("FITS Data Extension Value not in Parameter File")
    extn<-0
  }

  #Extension number of Error Map in FITS Header
  ID="ErrorExtn"
  extnerr<-as.numeric(params[ID,1])
  if (is.na(extnerr)) {
    warning("FITS Error Extension Value not in Parameter File")
    extnerr<-0
  }

  #Extension number of Mask Map in FITS Header
  ID="MaskExtn"
  extnmask<-as.numeric(params[ID,1])
  if (is.na(extnmask)) {
    warning("FITS Mask Extension Flag not in Parameter File")
    extnmask<-0
  }
  
  #Extension number of Data in FITS Header
  ID="WgtExtn"
  extnwgt<-as.numeric(params[ID,1])
  if (is.na(extnwgt)) {
    warning("FITS Weight Map Extension Value not in Parameter File")
    extnwgt<-0
  }

  #Do we want to force use of point sources?
  ID="PointSources"
  forcepointsources<-as.numeric(params[ID,1])
  if (is.na(forcepointsources)) {
    warning("Force Point Source Flag not in Parameter File")
    forcepointsources<-FALSE
  } else {
    forcepointsources<-(forcepointsources==1)
  }

  #Error Map scale factor
  ID="EFactor"
  Efactor<-as.numeric(params[ID,1])
  if (is.na(Efactor)) {
    warning("Error Map Scale Factor not in Parameter File")
    Efactor<-1
  }

  #Flux Correction (Scale) Factor
  ID="FluxCorr"
  fluxcorr<-as.numeric(params[ID,1])
  if (is.na(fluxcorr)) {
  warning("Flux Correction Factor not in Parameter File")
  fluxcorr<-1
  }

  cutrad<-NULL
  ra0<-NULL
  dec0<-NULL
  imfitsoutname<-NULL
  immfitsoutname<-NULL
  imefitsoutname<-NULL
  #Do we want to crop the input image(s)?
  ID="CropImage"
  cropimage<-params[ID,1]
  if (is.na(cropimage)) {
    warning("Crop Image Flag not in Parameter File")
    cropimage<-FALSE
  } else { cropimage<-(cropimage==1) }

  if (cropimage) {
    #What will the cropped image(s) be named
    ID="CropFitsName"
    imfitsoutname<-params[ID,1]
    if (is.na(imfitsoutname)) {
      warning("Cropped Images Filename not in Parameter File")
      imfitsoutname<-"croppedimage"
    }
    immfitsoutname<-paste(imfitsoutname,"_mask.fits",sep="")
    imefitsoutname<-paste(imfitsoutname,"_err.fits",sep="")
    imfitsoutname<-paste(imfitsoutname,".fits",sep="")
    ID="CropImRA0"
    ra0<-as.numeric(params[ID,1])
    if (is.na(ra0)) {
      warning("Cropped Images RA centre not in Parameter File")
    }
    ID="CropImDec0"
    dec0<-as.numeric(params[ID,1])
    if (is.na(dec0)) {
      warning("Cropped Images Dec Centre not in Parameter File")
    }
    ID="CropImRad"
    cutrad<-as.numeric(params[ID,1])
    if (is.na(cutrad)) {
      warning("Cropped Images cropping Radius not in Parameter File")
    }
  }

  #Confusion noise factor (in Janskys)
  ID="Confusion_Jy"
  conf<-as.numeric(params[ID,1])
  if (is.na(conf)) {
    warning("Confusion Noise Factor not in Parameter File")
    conf<-0
  }

  #Number of Processors available for computations
  ID="nProcessors"
  ncores<-as.numeric(params[ID,1])
  if (is.na(ncores)) {
    warning("Number of Processors Value not in Parameter File")
    ncores<-1
  }
  registerDoParallel(cores=ncores)

  #Is there an offset between the Input Catalogue angles and
  #N0E90 Angular Coordinates?
  ID="AngularOffset"
  angoffset<-params[ID,1]
  if (is.na(angoffset)) {
    warning("Angular Offset Flag not in Parameter File")
    angoffset<-FALSE
  } else { angoffset<-(angoffset==1) }

  #Is the map in Jy per Beam?
  ID="MapJyPerBeam"
  Jybm<-as.numeric(params[ID,1])
  if (is.na(Jybm)) {
    warning("Jansky Per Beam Flag not in Parameter File")
    Jybm<-0
  }

  #Do we want to perform higher precision integrations of
  #Apertures by resampling around the edges?
  ID="SmoothAper"
  resampleaperture<-as.numeric(params[ID,1])
  if (is.na(resampleaperture)) {
    warning("Smooth Aperture Flag  not in Parameter File")
    resampleaperture<-1
  }
  if (resampleaperture) {
    #What resolution do we want to upscale by?
    ID="ResamplingRes"
    upres<-as.numeric(params[ID,1])
    if (is.na(upres)) {
      warning("Resampling Resolution Factor not in Parameter File")
      upres<-2
    }
    #How many iterations of upscale do we want?
    ID="ResamplingIters"
    itersteps<-as.numeric(params[ID,1])
    if (is.na(itersteps)) {
      warning("Resampling Iterations Value not in Parameter File")
      itersteps<-10
    }
  } else {
    #If not - set defaults (#iters=0 performs no resampling)
    upres<-2
    itersteps<-0
  }

  #Number of PSF FWHM's that are added to the widths of the aperture stamps (which are, by default, ceiling(1.05*ApMajAxis) wide)
  #If we are not convolving with the PSF, PSF FWHM==0, and no buffer is added onto the default.
  ID="PSFConfidence"
  confidence<-as.numeric(params[ID,1])
  if (is.na(confidence)) {
    warning("PSFConfidence Value not in Parameter File")
    confidence<-0.95
  }  else if (confidence<=0) {
    warning("PSFConfidence Value is less than or equal to zero. Value should be strictly 0<confidence<1. Setting to 0.95")
    confidence<-0.95
  }

  #Size of the aperture stamp as a multiple of the aperture major axis
  ID="ApStampWidth"
  defbuff<-as.numeric(params[ID,1])
  if (is.na(defbuff)) {
    warning("ApStampWidth Value not in Parameter File")
    defbuff<-1.05
  }  else if (defbuff<=1) {
    warning("ApStampWidth Value is less than or equal to Unity. Value must be strictly > 1. Setting to 1.05")
    defbuff<-1.05
  }

  #Do we want to output the source mask only?
  ID="SourceMaskOnly"
  sourcemaskonly<-params[ID,1]
  if (is.na(sourcemaskonly)) {
    warning("Output Source Mask Only Flag not in Parameter File")
    sourcemaskonly<-FALSE
  } else { sourcemaskonly<-(sourcemaskonly==1) }

  if (!(sourcemaskonly)) {
    #Do we want to output the source mask at all?
    ID="WriteSourceMask"
    sourcemaskout<-params[ID,1]
    if (is.na(sourcemaskout)) {
      warning("Output Source Mask Flag not in Parameter File")
      sourcemaskout<-FALSE
      sourcemask   <-FALSE
    } else { 
      sourcemaskout<-(sourcemaskout==1) 
      sourcemask   <-(sourcemaskout==1) 
    }
  } else { 
    sourcemaskout<-TRUE 
    sourcemask   <-TRUE 
  }

  #Do we want to output the All Apertures Mask
  aafilename<-NULL
  ID="WriteAAMask"
  makeaamask<-params[ID,1]
  if (is.na(makeaamask)) {
    warning("Make All Apertures Mask Flag not in Parameter File")
    makeaamask<-FALSE
  } else { makeaamask<-(makeaamask==1) }

  if ( makeaamask ) {
    #Name of the output All Apertures file
    ID="AllApersFile"
    aafilename<-params[ID,1]
    if (is.na(aafilename)) {
      warning("All Apertures Mask Output Filename not in Parameter File")
      aafilename<-"AllApertures_Mask.fits"
    }
  }

  #Do we want to output the Convolved Apertures mask
  fafilename<-NULL
  ID="WriteFAMask"
  makefamask<-params[ID,1]
  if (is.na(makefamask)) {
    warning("Make Convolved Apertures Mask Flag not in Parameter File")
    makefamask<-FALSE
  } else { makefamask<-(makefamask==1) }

  if ( makefamask ) {
    #Name of the output Convolved Apertures Mask
    ID="ConvApersFile"
    fafilename<-params[ID,1]
    if (is.na(fafilename)) {
      warning("Convolved Apertures Mask Output Filename not in Parameter File")
      fafilename<-"AllConvolvedApertures_Mask.fits"
    }
  }

  #Do we want to output the Deblended Convolved Apertures mask
  dfafilename<-NULL
  ID="WriteDFAMask"
  makedfamask<-params[ID,1]
  if (is.na(makedfamask)) {
    warning("Make Deblended Convolved Apertures Mask Flag not in Parameter File")
    makedfamask<-FALSE
  } else { makedfamask<-(makedfamask==1) }

  if ( makedfamask ) {
    #Name of the output Convolved Apertures Mask
    ID="DeblConvApersFile"
    dfafilename<-params[ID,1]
    if (is.na(dfafilename)) {
      warning("Convolved Apertures Mask Output Filename not in Parameter File")
      dfafilename<-"AllDeblConvolvedApertures_Mask.fits"
    }
  }

  #Do we want to output the Residual image?
  residmap<-NULL
  ID="WriteResidMap"
  makeresidmap<-params[ID,1]
  if (is.na(makeresidmap)) {
    warning("Make Residual Map Flag not in Parameter File")
    makeresidmap<-TRUE
  } else { makeresidmap<-(makeresidmap==1) }

  if ( makeresidmap ) {
    #Name of the output Residual image
    ID="ResidImageFile"
    residmap<-params[ID,1]
    if (is.na(residmap)) {
      warning("Residual Map Filename not in Parameter File")
      residmap<-"ResidualImage.fits"
    }
  }

  #Do we want to output the Flux table?
  tableoutname<-NULL
  ID="WriteTable"
  writetab<-params[ID,1]
  if (is.na(writetab)) {
    warning("Write Table Flag not in Parameter File")
    writetab<-TRUE
  } else { writetab<-(writetab==1) }

  if ( writetab ) {
    #Name of output Flux Table
    ID="TableName"
    tableoutname<-params[ID,1]
    if (is.na(tableoutname)) {
      warning("Output Table Filename not in Parameter File")
      tableoutname<-"dfaResults.csv"
    }
  }

  #Do we want updates on how long things are taking?
  if (quiet) {
  showtime<-FALSE
  } else {
    ID="ShowTime"
    showtime<-params[ID,1]
    if (is.na(showtime)) {
      warning("ShowTime Flag not in Parameter File")
      showtime<-FALSE
    } else { showtime<-(showtime==1) }
  }

  #Do we want Diagnostic Output in Log File
  ID="Interactive"
  interact<-params[ID,1]
  if (is.na(interact)) {
    warning("Interactive Flag not in Parameter File")
    interact<-FALSE
  } else { interact<-(interact==1) }
  
  #What limit do we want for the use of masks what cross the mask edges
  ID="UseMaskLim"
  useMaskLim<-params[ID,1]
  if (is.na(useMaskLim)) {
    warning("UseMaskLim Value not in Parameter File")
    useMaskLim<-0.95
  }

  #Do we want Diagnostic Output in Log File
  ID="Diagnostic"
  diagnostic<-params[ID,1]
  if (is.na(diagnostic)) {
    warning("Diagnostic Flag not in Parameter File")
    diagnostic<-FALSE
  } else { diagnostic<-(diagnostic==1) }

  #Do we want Verbose Output in Log File
  ID="Verbose"
  verbose<-params[ID,1]
  if (is.na(verbose)) {
    warning("Verbose Flag not in Parameter File")
    verbose<-FALSE
  } else { verbose<-(verbose==1) }

  #Do we want a sample of the apertures to be output?
  ID="PlotSample"
  plotsample<-params[ID,1]
  if (is.na(plotsample)) {
    warning("Plot Sample Flag not in Parameter File")
    plotsample<-FALSE
  } else { plotsample<-(plotsample==1) }

  #Make Magnitudes in Output?
  ID="Magnitudes"
  Magnitudes<-params[ID,1]
  if (is.na(Magnitudes)) {
    warning("Make Magnitudes Flag not present in the Parameter File")
    Magnitudes<-TRUE
  } else { Magnitudes<-(Magnitudes==1) }

  #AB Vega Magnitude
  ID="ABVegaFlux"
  ABvegaflux<-as.numeric(params[ID,1])
  if (is.na(ABvegaflux)) {
    warning("AB Vega Flux Value not present in the Parameter File; using 1.0")
    ABvegaflux<-1.0
  }

  #Magnitudes Zero Point
  ID="MagZeroPoint"
  magZP<-as.numeric(params[ID,1])
  if (is.na(magZP)) {
    warning("Magnitudes Zero Point not present in the Parameter File; using 0.0")
    magZP<-0.0
  }

  #Magnitudes Zero Point
  ID="MagZPlabel"
  magZPlabel<-params[ID,1]
  if (is.na(magZPlabel)) {
    warning("FITS Zero Point Label not present in the Parameter File")
    magZPlabel<-"MagZP"
  }

  #Perform a Sky estimation & subtraction?
  ID="DoSkyEst"
  doskyest<-params[ID,1]
  if (is.na(doskyest)) {
    warning("Sky Estimate Flag not present in the Parameter File")
    doskyest<-FALSE
  } else { doskyest<-(doskyest==1) }

  #Calculate the Sky RMS? 
  ID="GetSkyRMS"
  getskyrms<-params[ID,1]
  if (is.na(getskyrms)) {
    warning("Get Sky RMS Flag not present in the Parameter File")
    getskyrms<-TRUE
  } else { getskyrms<-(getskyrms==1) }
  
  if (doskyest||getskyrms) { 
    #Sourcemask needed for SkyEstimate. If not TRUE, set to TRUE
    if (!sourcemask) {
      warning("Source Mask creation being forced because Sky Estimate Flag is TRUE")
      sourcemask<-TRUE
    }
    #Number of iterations used in sky estimation
    ID="SkyEstIters"
    skycutiters<-params[ID,1]
    if (is.na(skycutiters)) {
      warning("Sky Estimate Iterations Value not present in the Parameter File")
      skycutiters<-5
    }
    #Sigma level used in sky cut
    ID="SkyEstProbCut"
    skyprobcut<-as.numeric(params[ID,1])
    if (is.na(skyprobcut)) {
      warning("Sky Estimate Simga Cut Level not present in the Parameter File")
      skyprobcut<-3
    } 
    #Level of Corellated noise in image
    ID="SkyCorrelNoise"
    correl.noise<-as.numeric(params[ID,1])
    if (is.na(correl.noise)) {
      warning("Sky Correlated Noise Level not present in the Parameter File")
      correl.noise<-1
    } 
  } else {
    skycutiters<-0
    skyprobcut<-0
    correl.noise<-0
  }
  
  smfilename<-NULL
  if ( sourcemask ) {
    if (sourcemaskout) {
      #Name of SourceMask that is output
      ID="SourceMaskFile"
      smfilename<-params[ID,1]
      if (is.na(smfilename)) {
        warning("Source Mask filename not in Parameter File")
        smfilename<-"SourceMask.fits"
      }
    }
    #Name of SourceMask that is output
    ID="SourceMaskConfLim"
    smConfidenceLim<-as.numeric(params[ID,1])
    if (is.na(smConfidenceLim)) {
      warning("Source Mask Confidence Limit not in Parameter File")
      smConfidenceLim<-0.95
    }
  }
  
  #Name of Logfile to be output
  ID="LogFile"
  logfile<-params[ID,1]
  if (is.na(logfile)) {
    warning("Output Log Filename not present in the Parameter File")
    logfile<-"LAMDBAR_Log.txt"
  }

  #Send Message output to logfile
  sinkfile<-file(paste(pathout,logfile,sep=""),open="wt")
  sink(sinkfile, type="message")

  # Assign variables to LAMBDAR workspace
  assign("aafilename"       , aafilename       , envir = env) # A
  assign("ABvegaflux"       , ABvegaflux       , envir = env) #
  assign("angoffset"        , angoffset        , envir = env) #
  assign("beamarea_SOM_as"  , beamarea_SOM_as  , envir = env) # B
  assign("conf"             , conf             , envir = env) # C
  assign("confidence"       , confidence       , envir = env) #
  assign("cropimage"        , cropimage        , envir = env) #
  assign("correl.noise"     , correl.noise     , envir = env) #
  assign("catalogue"        , catalogue        , envir = env) #
  assign("cutrad"           , cutrad           , envir = env) #
  assign("contamlab"        , contamlab        , envir = env) #
  assign("catalab"          , catalab          , envir = env) #
  assign("datamap"          , datamap          , envir = env) # D
  assign("doskyest"         , doskyest         , envir = env) #
  assign("defbuff"          , defbuff          , envir = env) #
  assign("diagnostic"       , diagnostic       , envir = env) #
  assign("dec0"             , dec0             , envir = env) #
  assign("declab"           , declab           , envir = env) #
  assign("dfafilename"      , dfafilename      , envir = env) #
  assign("errormap"         , errormap         , envir = env) # E
  assign("extn"             , extn             , envir = env) #
  assign("extnerr"          , extnerr          , envir = env) #
  assign("extnmask"         , extnmask         , envir = env) #
  assign("extnwgt"          , extnwgt          , envir = env) #
  assign("Efactor"          , Efactor          , envir = env) #
  assign("fluxcorr"         , fluxcorr         , envir = env) # F
  assign("fafilename"       , fafilename       , envir = env) #
  assign("forcepointsources", forcepointsources, envir = env) #
  assign("filtcontam"       , filtcontam       , envir = env) #
  assign("fluxweightlab"    , fluxweightlab    , envir = env) #
  assign("gauss_fwhm_as"    , gauss_fwhm_as    , envir = env) # G
  assign("getskyrms"        , getskyrms        , envir = env) # 
  assign("itersteps"        , itersteps        , envir = env) # I
  assign("interact"         , interact         , envir = env) #
  assign("immfitsoutname"   , immfitsoutname   , envir = env) #
  assign("imefitsoutname"   , imefitsoutname   , envir = env) #
  assign("imfitsoutname"    , imfitsoutname    , envir = env) #
  assign("Jybm"             , Jybm             , envir = env) # J
  assign("makeresidmap"     , makeresidmap     , envir = env) # KLM
  assign("makedfamask"      , makedfamask      , envir = env) #
  assign("Magnitudes"       , Magnitudes       , envir = env) #
  assign("magZP"            , magZP            , envir = env) #
  assign("magZPlabel"       , magZPlabel       , envir = env) #
  assign("makefamask"       , makefamask       , envir = env) #
  assign("makeaamask"       , makeaamask       , envir = env) #
  assign("maskmap"          , maskmap          , envir = env) #
  assign("nopsf"            , nopsf            , envir = env) # N
  assign("nocontammap"      , nocontammap      , envir = env) #
  assign("ncores"           , ncores           , envir = env) #
  assign("overlay"          , overlay          , envir = env) # O
  assign("overlaymap"       , overlaymap       , envir = env) #
  assign("pathroot"         , pathroot         , envir = env) # P
  assign("pathout"          , pathout          , envir = env) #
  assign("params"           , params           , envir = env) #
  assign("plotsample"       , plotsample       , envir = env) #
  assign("psfmap"           , psfmap           , envir = env) #
  assign("resampleaperture" , resampleaperture , envir = env) # QR
  assign("ra0"              , ra0              , envir = env) #
  assign("ralab"            , ralab            , envir = env) #
  assign("residmap"         , residmap         , envir = env) #
  assign("sourcemask"       , sourcemask       , envir = env) # S
  assign("sourcemaskonly"   , sourcemaskonly   , envir = env) #
  assign("smConfidenceLim"  , smConfidenceLim  , envir = env) #
  assign("sinkfile"         , sinkfile         , envir = env) #
  assign("showtime"         , showtime         , envir = env) #
  assign("skycutiters"      , skycutiters      , envir = env) #
  assign("skyprobcut"       , skyprobcut       , envir = env) #
  assign("semimajlab"       , semimajlab       , envir = env) #
  assign("semiminlab"       , semiminlab       , envir = env) #
  assign("smfilename"       , smfilename       , envir = env) #
  assign("tableoutname"     , tableoutname     , envir = env) # T
  assign("thetalab"         , thetalab         , envir = env) #
  assign("upres"            , upres            , envir = env) # U
  assign("useMaskLim"       , useMaskLim       , envir = env)
  assign("verbose"          , verbose          , envir = env) # V
  assign("writetab"         , writetab         , envir = env) # W
                                                              # XYZ

  #Print any warnings
  if ((!quiet)&(!(is.null(warnings())))) {
  cat('\n')
  print(warnings())
  cat('                  ')
  }

  #-----Diagnostic------#
  if (diagnostic) {
    if (!(is.null(warnings()))) {
    sink(sinkfile, type="output")
    print(warnings())
    sink(type="output")
    }
    message("Variables set at end of Parameter File read")
    sink(file=sinkfile, type="output", append=TRUE)
    print(environment())
    #print(ls.str(envir=env))
    sink(file=NULL, type="output")
  }

  #Remove unneeded variables
  rm(params)

  #Finished Setup of Parameter Space
  if (!quiet) { cat(" - Done\n") }

}

readparfile <-
function(parfile=NA, starttime=NA, quiet=FALSE, env=NULL){
  #Procedure to setup Parameter Space (Read .par File) {{{
  if (is.null(env)) { env<-environment() }
  #}}}

  #Check Calling Syntax {{{
  if (is.na(parfile)) {
    stop("Parameter file not supplied. To create the default parameter file, run MeasureFluxes('--makepar').")
  }
  if (is.na(starttime)) {
    warning("Start time not supplied - using current clock time",call.=FALSE)
    starttime=proc.time()[3]
  }
  #}}}

  #Print Header {{{
  if (!quiet) { cat(paste('------------------------------------------------------\n'))
                cat(paste('   LAMBDAR     version ',packageVersion("LAMBDAR"),': MeasureFluxes\n'))
                cat(paste('------------------------------------------------------\n'))
                cat('Initialising Workspace {\n')
                cat('   Reading Parameter File   ') }
  #}}}

  #Test Reading of Parameter File {{{
  no.params<-try(max(count.fields(parfile)))
  if (class(no.params)=="try-error") {
  #Stop on Error
    if (grepl("cannot open the connection",no.params[1])) {
      cause="File not found"
    } else {
      cause="Cause unknown"
    }
    stop(paste("Parameter file read failed:",cause))
  }
  params<-try(read.table(parfile, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#", row.names=1, fill=TRUE, col.names=1:no.params), silent=TRUE)
  if (class(params)=="try-error") {
    #Stop on Error
    if (grepl("duplicate 'row.names'",params[1])) {
      cause="Duplicate Parameters in Parameter File"
    } else if (grepl("cannot open the connection",params[1])) {
      cause="File not found"
    } else {
      cause="Cause unknown"
    }
    stop(paste("Parameter file read failed:",cause))
  }
  if (!quiet) { cat('- Done\n   Assigning Parameter Variables') }
  #}}}

  ## Assign Parameter File values to variables {{{
  #Root Directory path {{{
  ID="RootDirectory"
  pathroot<-params[ID,1]
  if ((length(pathroot)==0)||(is.na(pathroot))) {
    stop("Root Directory Path not in Parameter File")
  }
  #Ensure path ends in a '/'
  if (lastnchar(pathroot,1) != '/') { pathroot<-paste(pathroot,'/',sep="") }
  #}}}

  #Root Directory path {{{
  ID="WorkingDirectory"
  ind<-which(params[ID,]!="")
  pathwork<-params[ID,ind]
  if ((length(ind)==0)||(is.na(pathroot))) {
    stop("Working Directory Path not in Parameter File")
  } else {
    pathwork<-try(suppressWarnings(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
    if (class(pathwork)=="try-error") {
      pathwork<-params[ID,ind]
    }
  }
  #Ensure path ends in a '/'
  if (any(lastnchar(pathwork,1) != '/')) { pathwork<-paste(pathwork,'/',sep="") }
  #}}}

  #Output Directory path {{{
  ID="OutputDirectory"
  ind<-which(params[ID,]!="")
  pathout<-params[ID,ind]
  if ((length(ind)==0)||(is.na(pathout))) {
    stop("Output Path not in Parameter File")
  } else {
    pathout<-try(suppressWarnings(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
    if (class(pathout)=="try-error") {
      pathout<-params[ID,ind]
    }
  }
  #Ensure path ends in a '/'
  if (lastnchar(pathout,1) != '/') { pathout<-paste(pathout,'/',sep="") }
  #}}}

  #Beam area in square arcsec {{{
  ID="BeamArea_SqAS"
  ind<-which(params[ID,]!="")
  beamarea_SOM_as<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(beamarea_SOM_as))) {
    warning("Beamarea Parameter not present in Parameter File",call.=FALSE)
    beamarea_SOM_as<-0
  }
  #}}}

  #Do we want to Convolve the apertures with a PSF {{{
  ID="PSFConvolve"
  ind<-which(params[ID,]!="")
  psffilt<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(psffilt))) {
    if ((length(ind)==1)) {
      psffilt<-try(as.numeric(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
      if (class(psffilt)=="try-error") {
        warning("PSFConvolve Parameter Table read failed. Using Default",call.=FALSE)
        psffilt<-0
      }
      if (is.na(psffilt)) {
        warning("PSFConvolve Parameter not in Parameter File. Using Default",call.=FALSE)
        psffilt<-0
      }
    } else {
      warning("PSF Convolve Flag not in Parameter File",call.=FALSE)
      psffilt<-0
    }
  }
  psffilt<-psffilt==1
  #}}}

  #Do we want PSF Matched Photometry? {{{
  ID="PSFMatched"
  ind<-which(params[ID,]!="")
  PSFMatched<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(PSFMatched))) {
    if ((length(ind)==1)) {
      PSFMatched<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(PSFMatched)=="try-error") {
        warning("Randoms Correction Flag not in Parameter File",call.=FALSE)
        PSFMatched<-0
      }
      if (is.na(PSFMatched)) {
        warning("Randoms Correction Flag not in Parameter File",call.=FALSE)
        PSFMatched<-0
      }
    } else {
      warning("Randoms Correction Flag not in Parameter File",call.=FALSE)
      PSFMatched<-0
    }
  }
  PSFMatched<-(PSFMatched==1)
  #}}}

  #PSF map filename {{{
  ID="PSFMap"
  ind<-which(params[ID,]!="")
  psfmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(psfmap))) {
    warning("PSF Map Filename not in Paramter File",call.=FALSE)
    psfmap<-"NONE"
  }
  #Determine if provided psfmap is an image or filelist {{{
  if ((length(psfmap)==1)&(psfmap!="NONE")&(!grepl(".fits", psfmap))) {
    #One file provided without .fits extension - must be filelist
    psfmap<-try(c(t(read.table(file.path(pathroot,psfmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(psfmap)=="try-error") {
      warning("PSF Map Filename not in Paramter File",call.=FALSE)
      psfmap<-"NONE"
    }
    if (is.na(psfmap)) {
      warning("PSF Map Filename not in Paramter File",call.=FALSE)
      psfmap<-"NONE"
    }
  }
  #}}}
  #}}}

  #If no PSF map, get gaussian FWHM {{{
  if (any(psfmap=="NONE")) {
    #FWHM of seeing gaussian
    ID="Gauss_FWHM_AS"
    ind<-which(params[ID,]!="")
    gauss_fwhm_as<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(gauss_fwhm_as))) {
      if ((length(ind)==1)) {
        gauss_fwhm_as<-try(as.numeric(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
        if (class(gauss_fwhm_as)=="try-error") {
          warning("Gauss_FWHM_AS Parameter Table read failed. Using Default",call.=FALSE)
          gauss_fwhm_as<-0.0
        }
        if (is.na(gauss_fwhm_as)) {
          warning("Gauss_FWHM_AS Parameter not in Parameter File. Using Default",call.=FALSE)
          gauss_fwhm_as<-0.0
        }
      } else {
        gauss_fwhm_as<-0.0
      }
    }
    #Make sure PSF maps and Gauss FWHM vals are conformable {{{
    if (length(gauss_fwhm_as)!=length(psfmap)) {
      gauss_fwhm_as<-rep(gauss_fwhm_as[1], length(psfmap))
    }#}}}

    #Make sure files with PSF maps have Gauss FWHM vals set to 0 {{{
    ind<-which(psfmap!="NONE")
    if (length(ind)>0) { gauss_fwhm_as[ind]<-0.0 }
    #}}}

    #If we want convolution, there is no PSF, and no gaussian FWHM provided - ERROR {{{
    ind<-which(psfmap=="NONE")
    if ((psffilt)&(any(gauss_fwhm_as[ind]==0.0))) {
      cat(" - Error\n")
      str<-paste("Loops with bad  parameters:",paste(which(psfmap=="NONE" & gauss_fwhm_as==0.0),collapse=", ",sep=""))
      stop(paste("Parameter file does not provide either PSF map or Gaussian FWHM for One or more files.\n",str,sep=""))
    }
    #}}}
  } else {
    #If there is a PSF Map - set gauss fwhm to zero {{{
    gauss_fwhm_as<-0.0
    #}}}
  }
  #}}}

  #Flag loops with no provided PSF {{{
  nopsf<-((psfmap=="NONE") & (gauss_fwhm_as==0.0))
  #}}}

  #Check for problematic inputs {{{
  if (any(PSFMatched & !psffilt)) {
    #message("WARNING: You've asked for PSF Matched Photometry, and yet specified No PSF Convolution... You could be in for some funky results",call.=FALSE)
    warning("You've asked for PSF Matched Photometry, and yet specified No PSF Convolution... You could be in for some funky results",call.=FALSE)
  }
  #}}}

  #Perform Contaminant removal? {{{
  nocontammap<-NULL
  ID="RemoveContam"
  ind<-which(params[ID,]!="")
  filtcontam<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(filtcontam))) {
    if ((length(ind)==1)) {
      filtcontam<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(filtcontam)=="try-error") {
        warning("Remove Contaminants Flag not in Parameter File",call.=FALSE)
        filtcontam<-FALSE
      } else { filtcontam<-(filtcontam==1) }
      if (is.na(filtcontam)) {
        warning("Remove Contaminants Flag not in Parameter File",call.=FALSE)
        filtcontam<-FALSE
      }
    } else {
      warning("Remove Contaminants Flag not in Parameter File",call.=FALSE)
      filtcontam<-FALSE
    }
  } else { filtcontam<-filtcontam==1 }
  #}}}

  #Contaminant Image Filename {{{
  if ( filtcontam ) {
    #Name of the output Residual image
    ID="NoContamImageFile"
    ind<-which(params[ID,]!="")
    nocontammap<-params[ID,ind]
    if ((length(ind)==0)||is.na(nocontammap)) {
      warning("Contaminant Subtracted Residual Map Filename not in Parameter File",call.=FALSE)
      nocontammap<-"NoContamResidualImage.fits"
    }
  }
  #}}}

  #Name of Source Catalogue {{{
  ID="Catalogue"
  ind<-which(params[ID,]!="")
  catalogue<-params[ID,ind]
  if ((length(ind)==0)||(is.na(catalogue))) {
    stop("Catalogue Path not in Parameter File")
  }
  #Determine if provided weightmap is an image or filelist {{{
  if ((length(catalogue)==1)&(!grepl(".csv",catalogue))&(!grepl(".Rdata",catalogue))&(!grepl(".fits", catalogue))) {
    #One file provided without .fits extension - must be filelist
    catalogue<-try(c(t(read.table(file.path(pathroot,catalogue), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(catalogue)=="try-error") {
      #Stop on Error
      stop("Catalogue Filelist read failed")
    }
    if (is.na(catalogue)) {
      stop("Catalogue Filelist read failed")
    }
  }
  #}}}
  #}}}

  #What is the title of the Catalogue's ID Column? {{{
  ID="CatIDColumnLabel"
  catalab<-params[ID,1]
  if (is.na(catalab)) {
    warning("Catalogue CATID Column Label not in Parameter File; using 'CATAID'",call.=FALSE)
    catalab<-"CATAID"
  }
  #}}}

  #What is the title of the Catalogue's RA Column? {{{
  ID="RAColumnLabel"
  ralab<-params[ID,1]
  if (is.na(ralab)) {
    warning("Catalogue RA Column Label not in Parameter File; using 'ALPHA_J2000'",call.=FALSE)
    ralab<-"ALPHA_J2000"
  }#}}}

  #What is the title of the Catalogue's Dec Column? {{{
  ID="DecColumnLabel"
  declab<-params[ID,1]
  if (is.na(declab)) {
    warning("Catalogue Dec Column Label not in Parameter File; using 'DELTA_J2000'",call.=FALSE)
    declab<-"DELTA_J2000"
  }#}}}

  #What is the title of the Catalogue's Theta Column? {{{
  ID="ThetaColumnLabel"
  thetalab<-params[ID,1]
  if (is.na(thetalab)) {
    warning("Catalogue Theta Column Label not in Parameter File; using 'THETA_J2000'",call.=FALSE)
    thetalab<-"THETA_J2000"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMajColumnLabel"
  semimajlab<-params[ID,1]
  if (is.na(semimajlab)) {
    warning("Catalogue SemiMajor Axis Column Label not in Parameter File; using 'SEMIMAJ_AS'",call.=FALSE)
    semimajlab<-"SEMIMAJ_AS"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMinColumnLabel"
  semiminlab<-params[ID,1]
  if (is.na(semiminlab)) {
    warning("Catalogue SemiMinor Axis Column Label not in Parameter File; using 'SEMIMIN_AS'",call.=FALSE)
    semiminlab<-"SEMIMIN_AS"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="ContamColumnLabel"
  contamlab<-params[ID,1]
  if (is.na(contamlab)) {
    warning("Catalogue SemiMajor Axis Column Label not in Parameter File; using 'CONTAM'",call.=FALSE)
    contamlab<-"CONTAM"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="FluxWgtColumnLabel"
  fluxweightlab<-params[ID,1]
  if (is.na(fluxweightlab)) {
    warning("Catalogue FluxWeight Column Label not in Parameter File; using 'FLUXWEIGHT'",call.=FALSE)
    fluxweightlab<-"FLUXWEIGHT"
  }#}}}

  #Name of Data Image {{{
  ID="DataMap"
  ind<-which(params[ID,]!="")
  datamap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(datamap))) {
    stop("Data Map Path not in Parameter File")
  }
  #Determine if provided datamap is an image or filelist {{{
  if ((length(datamap)==1)&(datamap!="NONE")&(!grepl(".fits", datamap))) {
    #One file provided without .fits extension - must be filelist
    datamap<-try(c(t(read.table(file.path(pathroot,datamap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(datamap)=="try-error") {
      #Stop on Error
      stop("Datamap Filelist read failed")
    }
    if (is.na(datamap)) {
      stop("Data Map Path not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Name of Error Map {{{
  ID="ErrorMap"
  ind<-which(params[ID,]!="")
  errormap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(errormap))) {
    warning("Error Map Path not in Parameter File",call.=FALSE)
    errormap<-"NONE"
  }
  #Determine if provided errormap is an image or filelist {{{
  if ((length(errormap)==1)&(is.na(as.numeric(errormap)))&(errormap!="NONE")&(!grepl(".fits", errormap))) {
    #One file provided without .fits extension - must be filelist
    errormap<-try(c(t(read.table(file.path(pathroot,errormap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(errormap)=="try-error") {
      #Stop on Error
      stop("ErrorMap Filelist read failed")
    }
    if (is.na(errormap)) {
      stop("ErrorMap Path not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Name of Mask Map {{{
  ID="MaskMap"
  ind<-which(params[ID,]!="")
  maskmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(maskmap))) {
    warning("Mask Map Path not in Parameter File",call.=FALSE)
    maskmap<-"NONE"
  }
  #Determine if provided maskmap is an image or filelist {{{
  if ((length(maskmap)==1)&(maskmap!="NONE")&(!grepl(".fits", maskmap))) {
    #One file provided without .fits extension - must be filelist
    maskmap<-try(c(t(read.table(file.path(pathroot,maskmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(maskmap)=="try-error") {
      #Stop on Error
      stop("maskmap Filelist read failed")
    }
    if (is.na(maskmap)) {
      stop("Data Map Path not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Name of Weight Map {{{
  ID="WeightMap"
  ind<-which(params[ID,]!="")
  wgtmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(wgtmap))) {
    warning("Weight Map Path not in Parameter File",call.=FALSE)
    wgtmap<-"NONE"
  }
  #Determine if provided weightmap is an image or filelist {{{
  if ((length(wgtmap)==1)&(wgtmap!="NONE")&(!grepl(".fits", wgtmap))) {
    #One file provided without .fits extension - must be filelist
    wgtmap<-try(c(t(read.table(file.path(pathroot,wgtmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(wgtmap)=="try-error") {
      #Stop on Error
      stop("weightmap Filelist read failed")
    }
    if (is.na(wgtmap)) {
      stop("Data Map Path not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Zero Point of Weight Map {{{
  ID="WeightMapZP"
  ind<-which(params[ID,]!="")
  wgtzp<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(wgtzp))) {
    warning("Weight Map Zero Point not in Parameter File; Using 0",call.=FALSE)
    wgtzp<-0
  }#}}}

  #Extension number of Data in FITS Header {{{
  ID="DataExtn"
  ind<-which(params[ID,]!="")
  extn<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extn))) {
    if ((length(ind)==1)) {
      extn<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(extn)=='try-error') {
        warning("DataExtn Parameter Table read failed. Using Default",call.=FALSE)
        extn<-0
      }
      if (is.na(extn)) {
        warning("DataExtn Parameter not in Parameter File. Using Default",call.=FALSE)
        extn<-0
      }
    } else {
      extn<-0
    }
  }
  #}}}

  #Extension number of Error Map in FITS Header {{{
  ID="ErrorExtn"
  ind<-which(params[ID,]!="")
  extnerr<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extnerr))) {
    if ((length(ind)==1)) {
      extnerr<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(extnerr)=='try-error') {
        warning("ErrorExtn Parameter Table read failed. Using Default",call.=FALSE)
        extnerr<-0
      }
      if (is.na(extnerr)) {
        warning("ErrorExtn Parameter not in Parameter File. Using Default",call.=FALSE)
        extnerr<-0
      }
    } else {
      extnerr<-0
    }
  }
  #}}}

  #Extension number of Mask Map in FITS Header {{{
  ID="MaskExtn"
  ind<-which(params[ID,]!="")
  extnmask<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extnmask))) {
    if ((length(ind)==1)) {
      extnmask<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(extnmask)=='try-error') {
        warning("MaskExtn Parameter Table read failed. Using Default",call.=FALSE)
        extnmask<-0
      }
      if (is.na(extnmask)) {
        warning("MaskExtn Parameter not in Parameter File. Using Default",call.=FALSE)
        extnmask<-0
      }
    } else {
      extnmask<-0
    }
  }
  #}}}

  #Extension number of Data in FITS Header {{{
  ID="WeightExtn"
  ind<-which(params[ID,]!="")
  extnwgt<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extnwgt))) {
    if ((length(ind)==1)) {
      extnwgt<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(extnwgt)=='try-error') {
        warning("WeightExtn Parameter Table read failed. Using Default",call.=FALSE)
        extnwgt<-0
      }
      if (is.na(extnwgt)) {
        warning("WeightExtn Parameter not in Parameter File. Using Default",call.=FALSE)
        extnwgt<-0
      }
    } else {
      extnwgt<-0
    }
  }
  #}}}

  #Do we want to force use of point sources? {{{
  ID="PointSources"
  ind<-which(params[ID,]!="")
  forcepointsources<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(forcepointsources)) {
    if ((length(ind)==1)) {
      forcepointsources<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(forcepointsources)=='try-error') {
        warning("PointSources Parameter Table read failed. Using Default",call.=FALSE)
        forcepointsources<-0
      }
      if (is.na(forcepointsources)) {
        warning("PointSources Parameter not in Parameter File. Using Default",call.=FALSE)
        forcepointsources<-0
      }
    } else {
      forcepointsources<-0
    }

  }
  forcepointsources<-(forcepointsources==1)
  #}}}

  #Check for silly inputs {{{
  if (any(PSFMatched & !forcepointsources)) {
    #message("WARNING: You've asked for PSF Matched Photometry, and not forced Point Sources to be used. PSF matched is designed for Point Sources only, so this could behave poorly.",call.=FALSE)
    warning("You've asked for PSF Matched Photometry, and not forced Point Sources to be used. PSF matched is designed for Point Sources only, so this could behave poorly.",call.=FALSE)
  }
  #}}}


  #Error Map scale factor #{{{
  ID="EFactor"
  ind<-which(params[ID,]!="")
  Efactor<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(Efactor))) {
    if ((length(ind)==1)) {
      Efactor<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(Efactor)=='try-error') {
        warning("EFactor Parameter Table read failed. Using Default",call.=FALSE)
        Efactor<-1
      }
      if (is.na(Efactor)) {
        warning("EFactor Parameter not in Parameter File. Using Default",call.=FALSE)
        Efactor<-1
      }
    } else {
      Efactor<-1
    }
  }
  #}}}

  #Flux Correction (Scale) Factor {{{
  ID="FluxCorr"
  ind<-which(params[ID,]!="")
  fluxcorr<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(fluxcorr))) {
    if ((length(ind)==1)) {
      fluxcorr<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(fluxcorr)=='try-error') {
        warning("FluxCorr Parameter Table read failed. Using Default",call.=FALSE)
        fluxcorr<-1
      }
      if (is.na(fluxcorr)) {
        warning("FluxCorr Parameter not in Parameter File. Using Default",call.=FALSE)
        fluxcorr<-1
      }
    } else {
      fluxcorr<-1
    }
  }
  #}}}

  #Initialise Cropping parameters {{{
  cutrad<-NULL
  ra0<-NULL
  dec0<-NULL
  imfitsoutname<-NULL
  immfitsoutname<-NULL
  imwgtfitsoutname<-NULL
  imefitsoutname<-NULL
  #}}}

  #Do we want to crop the input image(s)? {{{
  ID="CropImage"
  ind<-which(params[ID,]!="")
  cropimage<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(cropimage))) {
    if ((length(ind)==1)) {
      cropimage<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(cropimage)=='try-error') {
        warning("CropImage Parameter Table read failed. Using Default",call.=FALSE)
        cropimage<-0
      }
      if (is.na(cropimage)) {
        warning("CropImage Parameter not in Parameter File. Using Default",call.=FALSE)
        cropimage<-0
      }
    } else {
      warning("CropImage Flag not in Parameter File",call.=FALSE)
      cropimage<-0
    }
  }
  cropimage<-(cropimage==1)
  #}}}

  #Cropped image parameters {{{
  if (cropimage) {
    #What will the cropped image(s) be named {{{
    ID="CropFitsName"
    imfitsoutname<-params[ID,1]
    if (is.na(imfitsoutname)) {
      warning("Cropped Images Filename not in Parameter File",call.=FALSE)
      imfitsoutname<-"croppedimage"
    }
    immfitsoutname<-paste(imfitsoutname,"_mask.fits",sep="")
    imwgtfitsoutname<-paste(imfitsoutname,"_wgt.fits",sep="")
    imefitsoutname<-paste(imfitsoutname,"_err.fits",sep="")
    imfitsoutname<-paste(imfitsoutname,".fits",sep="")
    #}}}
    #Cropped image RA {{{
    ID="CropImRA0"
    ind<-which(params[ID,]!="")
    ra0<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(ra0))) {
      #Try Reading Table:
      ra0<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (is.na(ra0)||(class(ra0)=="try-error")) {
        #Warn on Error
        warning("Cropped Images RA centre not in Parameter File",call.=FALSE)
        ra0<- -999
      }
    }
    #}}}
    #Cropped image Dec {{{
    ID="CropImDec0"
    ind<-which(params[ID,]!="")
    dec0<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(dec0))) {
      #Try Reading Table:
      dec0<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (is.na(dec0)||(class(dec0)=="try-error")) {
        #Warn on Error
        warning("Cropped Images Dec Centre not in Parameter File",call.=FALSE)
        dec0<- -999
      }
    }
    #}}}
    #Cropped image radius {{{
    ID="CropImRad"
    ind<-which(params[ID,]!="")
    cutrad<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(cutrad))) {
      cutrad<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (is.na(cutrad)||(class(cutrad)=="try-error")) {
        warning("Cropped Images cropping Radius not in Parameter File",call.=FALSE)
        cutrad<-0.5
      }
    }
    #}}}
  }
  #}}}

  #Confusion noise factor (in Janskys) {{{
  ID="Confusion_Jy"
  ind<-which(params[ID,]!="")
  conf<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(conf))) {
    if ((length(ind)==1)) {
      conf<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(conf)=='try-error') {
        warning("Confusion_Jy Parameter Table read failed. Using Default.",call.=FALSE)
        conf<-1
      }
      if (is.na(conf)) {
        warning("Confusion_Jy Parameter not in Parameter File. Using Default.",call.=FALSE)
        conf<-1
      }
    } else {
      warning("Confusion_Jy Parameter not in Parameter File. Using Default.",call.=FALSE)
      conf<-1
    }
  }
  #}}}

  #Number of Processors available for computations {{{
  ID="nProcessors"
  ncores<-as.numeric(params[ID,1])
  if (is.na(ncores)) {
    warning("Number of Processors Value not in Parameter File",call.=FALSE)
    ncores<-1
  }
  #}}}

  #Angular Offset {{{
  #Is there an offset between the Input Catalogue angles and
  #N0E90 Angular Coordinates?
  ID="AngularOffset"
  ind<-which(params[ID,]!="")
  angoffset<-params[ID,ind]
  if ((length(ind)==0)||is.na(angoffset)) {
    if ((length(ind)==1)) {
      angoffset<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(angoffset)=='try-error') {
        warning("AngularOffset Parameter Table read failed. Using Default",call.=FALSE)
        angoffset<-0
      }
      if (is.na(angoffset)) {
        angoffset<-0
      }
    } else {
      angoffset<-0
    }
  }
  angoffset<-angoffset==1
  #}}}

  #Is the map in Jy per Beam? {{{
  ID="MapJyPerBeam"
  ind<-which(params[ID,]!="")
  Jybm<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(Jybm))) {
    if ((length(ind)==1)) {
      Jybm<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(Jybm)=='try-error') {
        warning("MapJyPerBeam Parameter Table read failed. Using Default",call.=FALSE)
        Jybm<-0
      }
      if (is.na(Jybm)) {
        warning("MapJyPerBeam Parameter not in Parameter File. Using Default",call.=FALSE)
        Jybm<-0
      }
    }
  }
  #}}}

  #Resample Apertures {{{
  #Do we want to perform higher precision integrations of
  #Apertures by resampling around the edges?
  ID="SmoothAper"
  ind<-which(params[ID,]!="")
  resampleaperture<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(resampleaperture))) {
    if ((length(ind)==1)) {
      resampleaperture<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(resampleaperture)=="try-error") {
        warning("SmoothAper Parameter Table read failed. Using Default",call.=FALSE)
        resampleaperture<-1
      }
      if (is.na(resampleaperture)) {
        warning("Aperture Resample Flag not in Parameter File",call.=FALSE)
        resampleaperture<-1
      }
    } else {
      warning("Aperture Resample Flag not in Parameter File",call.=FALSE)
      resampleaperture<-1
    }
  }
  resampleaperture<-resampleaperture==1
  #}}}

  #Resample Parameters {{{
  if (any(resampleaperture)) {
    #What resolution do we want to upscale by? {{{
    ID="ResamplingRes"
    ind<-which(params[ID,]!="")
    upres<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(upres))) {
      if ((length(ind)==1)) {
        upres<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(upres)=='try-error') {
          warning("ResamplingRes Parameter Table read failed. Using Default",call.=FALSE)
          upres<-3
        }
        if (is.na(upres)) {
          warning("ResamplingRes Parameter not in Parameter File. Using Default",call.=FALSE)
          upres<-3
        }
      } else {
        upres<-3
      }
    }
    #}}}
    #How many iterations of upscale do we want? {{{
    ID="ResamplingIters"
    ind<-which(params[ID,]!="")
    itersteps<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(itersteps))) {
      if ((length(ind)==1)) {
        itersteps<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(itersteps)=='try-error') {
          warning("ResamplingIters Parameter Table read failed. Using Default",call.=FALSE)
          itersteps<-5
        }
        if (is.na(itersteps)) {
          warning("ResamplingIters Parameter not in Parameter File. Using Default",call.=FALSE)
          itersteps<-5
        }
      } else {
        itersteps<-5
      }
    }
    #}}}
  } else {
    #If not - set defaults (#iters=0 performs no resampling) {{{
    upres<-2
    itersteps<-0
    #}}}
  }
  #}}}

  #PSF Confidence {{{
  #Number of PSF FWHM's that are added to the widths of the aperture stamps (which are, by default, ceiling(1.05*ApMajAxis) wide)
  #If we are not convolving with the PSF, PSF FWHM==0, and no buffer is added onto the default.
  ID="PSFConfidence"
  ind<-which(params[ID,]!="")
  confidence<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(confidence)) {
    if ((length(ind)==1)) {
      confidence<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(confidence)=='try-error') {
        warning("PSFConfidence Parameter Table read failed. Using Default",call.=FALSE)
        confidence<-1
      }
      if (is.na(confidence)) {
        warning("PSFConfidence Parameter not in Parameter File. Using Default",call.=FALSE)
        confidence<-1
      }
    } else {
      confidence<-1
    }
  }
  #}}}

  #Size of the aperture stamp as a multiple of the aperture major axis {{{
  ID="ApStampWidth"
  ind<-which(params[ID,]!="")
  defbuff<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(defbuff)) {
    if ((length(ind)==1)) {
      defbuff<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(defbuff)=='try-error') {
        warning("ApStampWidth Parameter Table read failed. Using Default",call.=FALSE)
        defbuff<-1.05
      }
      if (is.na(defbuff)) {
        warning("ApStampWidth Parameter not in Parameter File. Using Default",call.=FALSE)
        defbuff<-1.05
      }
    } else {
      defbuff<-1.05
    }
  }
  if (defbuff<1) {
    warning("ApStampWidth Value is less than or equal to Unity. Value must be strictly > 1. Setting to 1.05",call.=FALSE)
    defbuff<-1.05
  }
  #}}}

  #Do we want to output the source mask only? {{{
  ID="SourceMaskOnly"
  ind<-which(params[ID,]!="")
  sourcemaskonly<-as.numeric(params[ID,ind])
  if (length(ind)==0||is.na(sourcemaskonly)) {
    if (length(ind)==1) {
      sourcemaskonly<-try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
      if (class(sourcemaskonly)=="try-error") {
        warning("SourceMaskOnly Parameter not in Parameter File",call.=FALSE)
        sourcemaskonly<-0
      }
      if (is.na(sourcemaskonly)) {
        warning("SourceMaskOnly Parameter not in Parameter File",call.=FALSE)
        sourcemaskonly<-0
      }
    } else {
      warning("SourceMaskOnly Parameter not in Parameter File",call.=FALSE)
      sourcemaskonly<-0
    }
  }
  sourcemaskonly<-(sourcemaskonly==1)
  #}}}

  #Sourcemask filename {{{
  if (any(!sourcemaskonly)) {
    #Do we want to output the source mask at all?
    ID="WriteSourceMask"
    ind<-which(params[ID,]!="")
    sourcemaskout<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(sourcemaskout)) {
      if (length(ind)==1) {
        sourcemaskout<-try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
        if (class(sourcemaskout)=="try-error") {
          #Warn on Error
          warning("Output Source Mask Flag not in Parameter File",call.=FALSE)
          sourcemaskout<-0
        }
        if (is.na(sourcemaskout)) {
          warning("Output Source Mask Flag not in Parameter File",call.=FALSE)
          sourcemaskout<-0
        }
      } else {
        warning("Output Source Mask Flag not in Parameter File",call.=FALSE)
        sourcemaskout<-0
      }
    }
  } else {
    sourcemaskout<-TRUE
    sourcemask   <-TRUE
  }
  sourcemaskout<-(sourcemaskout==1)
  sourcemask   <-(sourcemaskout==1)
  #}}}

  #Do we want to output the All Apertures Mask {{{
  aafilename<-NULL
  ID="WriteAAMask"
  ind<-which(params[ID,]!="")
  makeaamask<-params[ID,ind]
  if ((length(ind)==0)||is.na(makeaamask)) {
    warning("Make All Apertures Mask Flag not in Parameter File",call.=FALSE)
    makeaamask<-FALSE
  } else { makeaamask<-(makeaamask==1) }
  #}}}

  #All Apertures mask file name {{{
  if ( makeaamask ) {
    #Name of the output All Apertures file
    ID="AllApersFile"
    ind<-which(params[ID,]!="")
    aafilename<-params[ID,ind]
    if ((length(ind)==0)||is.na(aafilename)) {
      warning("All Apertures Mask Output Filename not in Parameter File",call.=FALSE)
      aafilename<-"AllApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Convolved Apertures mask {{{
  fafilename<-NULL
  ID="WriteFAMask"
  ind<-which(params[ID,]!="")
  makefamask<-params[ID,ind]
  if ((length(ind)==0)||is.na(makefamask)) {
    warning("Make Convolved Apertures Mask Flag not in Parameter File",call.=FALSE)
    makefamask<-FALSE
  } else { makefamask<-(makefamask==1) }
  #}}}

  #Convolved Apertures Filename {{{
  if ( makefamask ) {
    #Name of the output Convolved Apertures Mask
    ID="ConvApersFile"
    ind<-which(params[ID,]!="")
    fafilename<-params[ID,ind]
    if ((length(ind)==0)||is.na(fafilename)) {
      warning("Convolved Apertures Mask Output Filename not in Parameter File",call.=FALSE)
      fafilename<-"AllConvolvedApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Deblended Convolved Apertures mask {{{
  dfafilename<-NULL
  ID="WriteDFAMask"
  ind<-which(params[ID,]!="")
  makedfamask<-params[ID,ind]
  if ((length(ind)==0)||is.na(makedfamask)) {
    warning("Make Deblended Convolved Apertures Mask Flag not in Parameter File",call.=FALSE)
    makedfamask<-FALSE
  } else { makedfamask<-(makedfamask==1) }
  #}}}

  #Deblended Convolved Apertures Mask filename {{{
  if ( makedfamask ) {
    #Name of the output Convolved Apertures Mask
    ID="DeblConvApersFile"
    ind<-which(params[ID,]!="")
    dfafilename<-params[ID,ind]
    if ((length(ind)==0)||is.na(dfafilename)) {
      warning("Convolved Apertures Mask Output Filename not in Parameter File",call.=FALSE)
      dfafilename<-"AllDeblConvolvedApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Residual image? {{{
  residmap<-NULL
  ID="WriteResidMap"
  ind<-which(params[ID,]!="")
  makeresidmap<-params[ID,ind]
  if ((length(ind)==0)||is.na(makeresidmap)) {
    warning("Make Residual Map Flag not in Parameter File",call.=FALSE)
    makeresidmap<-TRUE
  } else { makeresidmap<-(makeresidmap==1) }
  #}}}

  #Residual Image filename {{{
  if ( makeresidmap ) {
    #Name of the output Residual image
    ID="ResidImageFile"
    ind<-which(params[ID,]!="")
    residmap<-params[ID,ind]
    if ((length(ind)==0)||is.na(residmap)) {
      warning("Residual Map Filename not in Parameter File",call.=FALSE)
      residmap<-"ResidualImage.fits"
    }
  }
  #}}}

  #Do we want to output the Flux table? {{{
  tableoutname<-NULL
  ID="WriteTable"
  ind<-which(params[ID,]!="")
  writetab<-params[ID,ind]
  if ((length(ind)==0)||is.na(writetab)) {
    warning("Write Table Flag not in Parameter File",call.=FALSE)
    writetab<-TRUE
  } else { writetab<-(writetab==1) }
  #}}}

  #Table Filename {{{
  if ( writetab ) {
    #Name of output Flux Table
    ID="TableName"
    ind<-which(params[ID,]!="")
    tableoutname<-params[ID,ind]
    if ((length(ind)==0)||(is.na(tableoutname))) {
      #Warn on Error
      warning("Output Table Filename not in Parameter File",call.=FALSE)
      tableoutname<-"dfaResults"
    } else {
      if (length(ind)==1) {
        tableoutname<-try(suppressWarnings(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
        if (class(tableoutname)=="try-error") {
          #Warn on Error
          tableoutname<-params[ID,1]
        }
        if (is.na(tableoutname)) {
          tableoutname<-params[ID,1]
        }
      }
    }
  }
  #}}}

  #Do we want updates on how long things are taking? {{{
  if (quiet) {
    showtime<-FALSE
  } else {
    ID="ShowTime"
    ind<-which(params[ID,]!="")
    showtime<-params[ID,ind]
    if ((length(ind)==0)||is.na(showtime)) {
      warning("ShowTime Flag not in Parameter File",call.=FALSE)
      showtime<-FALSE
    } else { showtime<-(showtime==1) }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Interactive"
  ind<-which(params[ID,]!="")
  interact<-params[ID,ind]
  if ((length(ind)==0)||is.na(interact)) {
    warning("Interactive Flag not in Parameter File",call.=FALSE)
    interact<-FALSE
  } else { interact<-(interact==1) }
  #}}}

  #What limit do we want for the use of masks what cross the mask edges {{{
  ID="UseMaskLim"
  ind<-which(params[ID,]!="")
  useMaskLim<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(useMaskLim)) {
    if ((length(ind)==1)) {
      useMaskLim<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(useMaskLim)=='try-error') {
        warning("UseMaskLim Parameter Table read failed. Using Default",call.=FALSE)
        useMaskLim<-0.2
      }
      if (is.na(useMaskLim)) {
        warning("UseMaskLim Parameter not in Parameter File. Using Default",call.=FALSE)
        useMaskLim<-0.2
      }
    } else {
      useMaskLim<-0.2
    }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Diagnostic"
  ind<-which(params[ID,]!="")
  diagnostic<-params[ID,ind]
  if ((length(ind)==0)||is.na(diagnostic)) {
    warning("Diagnostic Flag not in Parameter File",call.=FALSE)
    diagnostic<-FALSE
  } else { diagnostic<-(diagnostic==1) }
  #}}}

  #Do we want Verbose Output in Log File {{{
  if (diagnostic) {
    verbose<-TRUE
  } else {
    ID="Verbose"
    ind<-which(params[ID,]!="")
    verbose<-params[ID,ind]
    if ((length(ind)==0)||is.na(verbose)) {
      warning("Verbose Flag not in Parameter File",call.=FALSE)
      verbose<-FALSE
    } else { verbose<-(verbose==1) }
  }
  #}}}

  #Do we want a sample of the apertures to be output? {{{
  ID="PlotSample"
  ind<-which(params[ID,]!="")
  plotsample<-params[ID,ind]
  if ((length(ind)==0)||is.na(plotsample)) {
    warning("Plot Sample Flag not in Parameter File",call.=FALSE)
    plotsample<-FALSE
  } else { plotsample<-(plotsample==1) }
  ID="PlotAll"
  ind<-which(params[ID,]!="")
  plotall<-params[ID,ind]
  if ((length(ind)==0)||is.na(plotall)) {
    warning("Plot All Flag not in Parameter File",call.=FALSE)
    plotall<-FALSE
  } else { plotall<-(plotall==1) }
  if (plotall) { plotsample<-TRUE }
  #}}}

  #Make Magnitudes in Output? {{{
  ID="Magnitudes"
  ind<-which(params[ID,]!="")
  Magnitudes<-params[ID,ind]
  if ((length(ind)==0)||is.na(Magnitudes)) {
    warning("Make Magnitudes Flag not present in the Parameter File",call.=FALSE)
    Magnitudes<-TRUE
  } else { Magnitudes<-(Magnitudes==1) }
  #}}}

  #Magnitude Details {{{
  if (Magnitudes) {
    #AB Vega Magnitude {{{
    ID="ABVegaFlux"
    ind<-which(params[ID,]!="")
    ABvegaflux<-as.numeric(params[ID,ind])
    if (is.na(ABvegaflux)) {
      warning("AB Vega Flux Value not present in the Parameter File; using 1.0",call.=FALSE)
      ABvegaflux<-1.0
    }
    #}}}

    #Magnitudes Zero Point {{{
    ID="MagZeroPoint"
    ind<-which(params[ID,]!="")
    magZP<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(magZP))) {
      if (length(ind)==1) {
        magZP<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(magZP)=="try-error") {
          #Warn on Error
          warning("Magnitudes Zero Point not present/bad in the Parameter File; using 0.0",call.=FALSE)
          magZP<-0.0
        }
        if (is.na(magZP)) {
          #Warn on Error
          warning("Magnitudes Zero Point not present/bad in the Parameter File; using 0.0",call.=FALSE)
          magZP<-0.0
        }
      } else {
        #Warn on Error
        warning("Magnitudes Zero Point not present/bad in the Parameter File; using 0.0",call.=FALSE)
        magZP<-0.0
      }
    }
    #}}}

    #Magnitudes Zero Point {{{
    ID="MagZPLabel"
    ind<-which(params[ID,]!="")
    magZPlabel<-params[ID,ind]
    if ((length(ind)==0)||(is.na(magZPlabel))) {
      warning("FITS Zero Point Label not present in the Parameter File",call.=FALSE)
      magZPlabel<-"MagZP"
    }
    #}}}
  } else {
    magZP<-NA
    magZPlabel<-"MagZP"
    ABvegaflux<-1.0
  }
  #}}}

  #Randoms Correction {{{
  #Do we want to perform a randoms Correction to the errors of object fluxes?
  ID="RanCor"
  ind<-which(params[ID,]!="")
  RanCor<-params[ID,ind]
  if (length(ind)!=0){ if(RanCor!="execute") { RanCor<-as.numeric(RanCor) } }
  if ((length(ind)==0)||(is.na(RanCor))) {
    if ((length(ind)==1)) {
      RanCor<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(RanCor)=="try-error") {
        warning("Randoms Correction Flag not in Parameter File",call.=FALSE)
        RanCor<-0
      }
      if (is.na(RanCor)) {
        warning("Randoms Correction Flag not in Parameter File",call.=FALSE)
        RanCor<-0
      }
    } else {
      warning("Randoms Correction Flag not in Parameter File",call.=FALSE)
      RanCor<-0
    }
  }
  RanCor<-RanCor==1
  #Number of Randoms per Object {{{
  ID="nRandoms"
  ind<-which(params[ID,]!="")
  nRandoms<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(nRandoms))) {
    if ((length(ind)==1)) {
      nRandoms<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(nRandoms)=="try-error") {
        warning("Number of Randoms Flag not in Parameter File",call.=FALSE)
        nRandoms<-10
      }
      if (is.na(nRandoms)) {
        warning("Number of Randoms Flag not in Parameter File",call.=FALSE)
        nRandoms<-10
      }
    } else {
      warning("Number of Randoms Flag not in Parameter File",call.=FALSE)
      nRandoms<-10
    }
  }
  #}}}
  #}}}

  #Perform a Sky estimation & subtraction? {{{
  ID="DoSkyEst"
  ind<-which(params[ID,]!="")
  doskyest<-as.numeric(params[ID,ind])
  if ((length(ind)==0) || is.na(doskyest)) {
    if (length(ind)==1) {
      doskyest<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(doskyest)=="try-error") {
        #Warn on Error
        warning("Sky Estimate Flag not present in the Parameter File",call.=FALSE)
        doskyest<-0
      }
      if (is.na(doskyest)) {
        #Warn on Error
        warning("Sky Estimate Flag not present in the Parameter File",call.=FALSE)
        doskyest<-0
      }
    } else {
        #Warn on Error
        warning("Sky Estimate Flag not present in the Parameter File",call.=FALSE)
        doskyest<-0
    }
  }
  doskyest<-(doskyest==1)
  #}}}

  #Calculate the Sky RMS? {{{
  ID="GetSkyRMS"
  ind<-which(params[ID,]!="")
  getskyrms<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(getskyrms)) {
    if (length(ind)==1) {
      getskyrms<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(getskyrms)=="try-error") {
        #Warn on Error
        warning("GetSkyRMS Flag not present in the Parameter File",call.=FALSE)
        getskyrms<-0
      }
      if (is.na(getskyrms)) {
        #Warn on Error
        warning("GetSkyRMS Flag not present in the Parameter File",call.=FALSE)
        getskyrms<-0
      }
    } else {
        #Warn on Error
        warning("GetSkyRMS Flag not present in the Parameter File",call.=FALSE)
        getskyrms<-0
    }
  }
  getskyrms<-getskyrms==1
  #}}}

  #Sky Estimate Paramters {{{
  if (any(doskyest|getskyrms)) {
    #Sourcemask needed for SkyEstimate. If not TRUE, set to TRUE {{{
    if (any(!sourcemask & !(doskyest|getskyrms))) {
      message("Source Mask creation being forced for all runs with Sky Estimate Flag TRUE")
      sourcemask<-doskyest|getskyrms|sourcemask
    }
    #}}}

    #Number of iterations used in sky estimation {{{
    ID="SkyEstIters"
    ind<-which(params[ID,]!="")
    skycutiters<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||is.na(skycutiters)) {
      warning("Sky Estimate Iterations Value not present in the Parameter File",call.=FALSE)
      skycutiters<-5
    }
    #}}}

    #Sigma level used in sky cut {{{
    ID="SkyEstProbCut"
    ind<-which(params[ID,]!="")
    skyprobcut<-as.numeric(params[ID,ind])
    if ((length(ind)==0)|| (is.na(skyprobcut))) {
      if ((length(ind)==1)) {
        skyprobcut<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(skyprobcut)=='try-error') {
          warning("SkyEstProbCut Parameter Table read failed. Using Default",call.=FALSE)
          skyprobcut<-3
        }
        if (is.na(skyprobcut)) {
          warning("SkyEstProbCut Parameter not in Parameter File. Using Default",call.=FALSE)
          skyprobcut<-3
        }
      } else {
        skyprobcut<-3
      }
    }
    #}}}

    #Default Sky Value if estimation fails {{{
    ID="SkyDefault"
    ind<-which(params[ID,]!="")
    skydefault<-params[ID,ind]
    if ((length(ind)==0)||is.na(skydefault)) {
      warning("Sky Default Value not present in the Parameter File; using median",call.=FALSE)
      skydefault<-"median"
    }
    #}}}

    #Level of Corellated noise in image {{{
    ID="SkyCorrelNoise"
    ind<-which(params[ID,]!="")
    correl.noise<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(correl.noise))) {
      if ((length(ind)==1)) {
        correl.noise<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(correl.noise)=='try-error') {
          warning("SkyCorrelNoise Parameter Table read failed. Using Default",call.=FALSE)
          correl.noise<-1
        }
        if (is.na(correl.noise)) {
          warning("SkyCorrelNoise Parameter not in Parameter File. Using Default",call.=FALSE)
          correl.noise<-1
        }
      } else {
        correl.noise<-1
      }
    }
    #}}}
  } else {
    skycutiters<-0
    skydefault<-0
    skyprobcut<-0
    correl.noise<-1
  }
  #}}}

  #Sourcemask parameter {{{
  smfilename<-NULL
  smConfidenceLim<-NULL
  if ( sourcemask ) {
    if (sourcemaskout) {
      #Name of SourceMask that is output
      ID="SourceMaskFile"
      ind<-which(params[ID,]!="")
      smfilename<-params[ID,ind]
      if ((length(ind)==0)|| (is.na(smfilename))) {
        warning("Source Mask filename not in Parameter File",call.=FALSE)
        smfilename<-"SourceMask.fits"
      }
    }
    #Name of SourceMask that is output
    ID="TransmissionMap"
    ind<-which(params[ID,]!="")
    TransmissionMap<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(TransmissionMap)) {
      if (length(ind)==1) {
        TransmissionMap<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(TransmissionMap)=="try-error") {
          #Warn on Error
          warning("Transmission Map Flag not in Parameter File",call.=FALSE)
          TransmissionMap<-0
        }
        if (is.na(TransmissionMap)) {
          #Warn on Error
          warning("Transmission Map Flag not in Parameter File",call.=FALSE)
          TransmissionMap<-0
        }
      } else {
        #Warn on Error
        warning("Transmission Map Flag not in Parameter File",call.=FALSE)
        TransmissionMap<-0
      }
    }
    TransmissionMap<-(TransmissionMap==1)
    #SourceMask Confidence Limit
    ID="SourceMaskConfLim"
    ind<-which(params[ID,]!="")
    smConfidenceLim<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(smConfidenceLim)) {
      if (length(ind)==1) {
        smConfidenceLim<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(smConfidenceLim)=="try-error") {
          #Warn on Error
          warning("Sourcemask Confidence Specification not in Parameter File",call.=FALSE)
          smConfidenceLim<-0.95
        }
        if (is.na(smConfidenceLim)) {
          #Warn on Error
          warning("Sourcemask Confidence Specification not in Parameter File",call.=FALSE)
          smConfidenceLim<-0.95
        }
      } else {
        #Warn on Error
        warning("Sourcemask Confidence Specification not in Parameter File",call.=FALSE)
        smConfidenceLim<-0.95
      }
    }
  } else {
    TransmissionMap=FALSE
    smConfidenceLim=NA
  }
  #}}}

  #Set Minimum Aperture Radius {{{
  ID="MinApRad"
  ind<-which(params[ID,]!="")
  MinApRad<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(MinApRad))) {
    if (length(ind)==1) {
      MinApRad<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(MinApRad)=="try-error") {
        #Warn on Error
        warning("Minimum Aperture Specification not in Parameter File",call.=FALSE)
        MinApRad<-0
      }
      if (is.na(MinApRad)) {
        #Warn on Error
        warning("Minimum Aperture Specification not in Parameter File",call.=FALSE)
        MinApRad<-0
      }
    } else {
      #Warn on Error
      warning("Minimum Aperture Specification not in Parameter File",call.=FALSE)
      MinApRad<-0
    }
  }
  #}}}

  #Do we want a Memory-Safe run? {{{
  ID="MemorySafe"
  ind<-which(params[ID,]!="")
  memSafe<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(memSafe)) {
    if ((length(ind)==1)) {
      memSafe<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(memSafe)=='try-error') {
        warning("MemorySafe Parameter Table read failed. Using Default",call.=FALSE)
        memSafe<-0
      }
      if (is.na(memSafe)) {
        memSafe<-0
        warning("MemorySafe Parameter not in Parameter File. Using Default",call.=FALSE)
      }
    } else {
      memSafe<-0
      warning("MemorySafe Parameter not in Parameter File. Using Default",call.=FALSE)
    }
  }
  #}}}

  #What limit value do we want for the Aperture generation? #{{{
  ID="ApertureConfLimit"
  ind<-which(params[ID,]!="")
  apLimit<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(apLimit))) {
    if (length(ind)==1) {
      apLimit<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(apLimit)=="try-error") {
        #Warn on Error
        warning("Aperture Confidence Limit Value not in Parameter File",call.=FALSE)
        apLimit<-0.9
      }
      if (is.na(apLimit)) {
        #Warn on Error
        warning("Aperture Confidence Limit Value not in Parameter File",call.=FALSE)
        apLimit<-0.9
      }
    } else {
      #Warn on Error
      warning("Aperture Confidence Limit Value not in Parameter File",call.=FALSE)
      apLimit<-0.9
    }
  }
  #}}}

  #Do we want to iterate the fluxes? {{{
  ID="IterateFluxes"
  ind<-which(params[ID,]!="")
  iterateFluxes<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(iterateFluxes)) {
    if (length(ind)==1) {
      iterateFluxes<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(iterateFluxes)=="try-error") {
        #Warn on Error
        warning("IterateFluxes Flag Table read failed. Using Default.",call.=FALSE)
        iterateFluxes<-0
      }
      if (is.na(iterateFluxes)) {
        #Warn on Error
        warning("IterateFluxes Flag not present in the Parameter File",call.=FALSE)
        iterateFluxes<-0
      }
    } else {
        #Warn on Error
        warning("IterateFluxes Flag not present in the Parameter File",call.=FALSE)
        iterateFluxes<-0
    }
  }
  iterateFluxes<-(iterateFluxes==1)
  #}}}

  #If Iterating, how many iterations do we want? {{{
  ID="nIterations"
  ind<-which(params[ID,]!="")
  nIterations<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(nIterations)) {
    if ((length(ind)==1)) {
      nIterations<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(nIterations)=='try-error') {
        warning("nIterations Parameter Table read failed. Using Default",call.=FALSE)
        nIterations<-2
      }
      if (is.na(nIterations)) {
        warning("nIterations Parameter not in Parameter File. Using Default",call.=FALSE)
        nIterations<-2
      }
    } else {
      warning("nIterations Parameter not in Parameter File. Using Default",call.=FALSE)
      nIterations<-2
    }
  }
  #}}}

  #What format are the input fluxweights? {{{
  ID="FluxWgtType"
  ind<-which(params[ID,]!="")
  weightType<-tolower(params[ID,ind])
  if ((length(ind)==0)||is.na(weightType)) {
    warning("Fluxweight Type Parameter not present in the Parameter File",call.=FALSE)
    weightType<-"scale"
  } else if (weightType!="flux" & weightType!="mag" & weightType!="scale") {
    stop("Fluxweight Type Parameter has unknown value; valid values are 'flux', 'mag', and 'scale'")
  }
  #}}}

  #Do we want to fluxweight using Pixel Fluxes? {{{
  ID="UsePixelFluxWgts"
  ind<-which(params[ID,]!="")
  usePixelFluxWeights<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(usePixelFluxWeights)) {
    if (length(ind)==1) {
      usePixelFluxWeights<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(usePixelFluxWeights)=="try-error") {
        #Warn on Error
        warning("UsePixelFluxWgts Table read failed. Using Default.",call.=FALSE)
        usePixelFluxWeights<-0
      }
      if (is.na(usePixelFluxWeights)) {
        #Warn on Error
        warning("UsePixelFluxWgts Flag not present in the Parameter File",call.=FALSE)
        usePixelFluxWeights<-0
      }
    } else {
        #Warn on Error
        warning("UsePixelFluxWgts Flag not present in the Parameter File",call.=FALSE)
        usePixelFluxWeights<-0
    }
  }
  usePixelFluxWeights<-(usePixelFluxWeights==1)
  if (usePixelFluxWeights) { weightType<-"scale" }
  #}}}

  #Name of Logfile to be output {{{
  ID="LogFile"
  logfile<-params[ID,1]
  if (is.na(logfile)) {
    warning("Output Log Filename not present in the Parameter File",call.=FALSE)
    logfile<-"LAMDBAR_Log.txt"
  }
  #}}}
  #}}}

  # Print any warnings {{{
  if (!is.null(warnings())) {
    cat(" {\nWarnings in Parameter File read:\n")
    print(warnings())
    cat("} ")
  }

  #}}}

  # Assign variables to LAMBDAR workspace {{{
  assign("aafilename"       , aafilename       , envir = env) # A
  assign("ABvegaflux"       , ABvegaflux       , envir = env) #
  assign("angoffset"        , angoffset        , envir = env) #
  assign("apLimit"          , apLimit          , envir = env) #
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
  assign("iterateFluxes"    , iterateFluxes    , envir = env) # I
  assign("nIterations"      , nIterations      , envir = env) # I
  assign("interact"         , interact         , envir = env) #
  assign("immfitsoutname"   , immfitsoutname   , envir = env) #
  assign("imwgtfitsoutname" , imwgtfitsoutname , envir = env) #
  assign("imefitsoutname"   , imefitsoutname   , envir = env) #
  assign("imfitsoutname"    , imfitsoutname    , envir = env) #
  assign("Jybm"             , Jybm             , envir = env) # J
  assign("logfile"          , logfile          , envir = env) # KL
  assign("makeresidmap"     , makeresidmap     , envir = env) # M
  assign("makedfamask"      , makedfamask      , envir = env) #
  assign("Magnitudes"       , Magnitudes       , envir = env) #
  assign("magZP"            , magZP            , envir = env) #
  assign("magZPlabel"       , magZPlabel       , envir = env) #
  assign("makefamask"       , makefamask       , envir = env) #
  assign("makeaamask"       , makeaamask       , envir = env) #
  assign("maskmap"          , maskmap          , envir = env) #
  assign("MinApRad"         , MinApRad         , envir = env) #
  assign("memSafe"          , memSafe          , envir = env) #
  assign("nopsf"            , nopsf            , envir = env) # N
  assign("nocontammap"      , nocontammap      , envir = env) #
  assign("ncores"           , ncores           , envir = env) #
  assign("nRandoms"         , nRandoms         , envir = env) #
  assign("pathroot"         , pathroot         , envir = env) # P
  assign("pathwork"         , pathwork         , envir = env) #
  assign("pathout"          , pathout          , envir = env) #
  assign("plotsample"       , plotsample       , envir = env) #
  assign("plotall"          , plotall          , envir = env) #
  assign("psfmap"           , psfmap           , envir = env) #
  assign("PSFMatched"       , PSFMatched       , envir = env) #
  assign("psffilt"          , psffilt          , envir = env) #
  assign("resampleaperture" , resampleaperture , envir = env) # QR
  assign("ra0"              , ra0              , envir = env) #
  assign("ralab"            , ralab            , envir = env) #
  assign("RanCor"           , RanCor           , envir = env) #
  assign("residmap"         , residmap         , envir = env) #
  assign("sourcemask"       , sourcemask       , envir = env) # S
  assign("sourcemaskout"    , sourcemaskout    , envir = env) #
  assign("sourcemaskonly"   , sourcemaskonly   , envir = env) #
  assign("smConfidenceLim"  , smConfidenceLim  , envir = env) #
  assign("showtime"         , showtime         , envir = env) #
  assign("skycutiters"      , skycutiters      , envir = env) #
  assign("skydefault"       , skydefault       , envir = env) #
  assign("skyprobcut"       , skyprobcut       , envir = env) #
  assign("semimajlab"       , semimajlab       , envir = env) #
  assign("semiminlab"       , semiminlab       , envir = env) #
  assign("smfilename"       , smfilename       , envir = env) #
  assign("tableoutname"     , tableoutname     , envir = env) # T
  assign("thetalab"         , thetalab         , envir = env) #
  assign("TransmissionMap"  , TransmissionMap  , envir = env) #
  assign("upres"            , upres            , envir = env) # U
  assign("useMaskLim"       , useMaskLim       , envir = env)
  assign("usePixelFluxWeights", usePixelFluxWeights , envir = env)
  assign("verbose"          , verbose          , envir = env) # V
  assign("writetab"         , writetab         , envir = env) # W
  assign("weightType"       , weightType       , envir = env) #
  assign("wgtmap"           , wgtmap           , envir = env) #
  assign("wgtzp"            , wgtzp            , envir = env) #
                                                              # XYZ
  #}}}

  #Remove unneeded variables {{{
  rm(params)
  #}}}

  #Finished Setup of Parameter Space {{{
  if (!quiet) { cat(" - Done\n") }
  return=NULL
  #}}}

}

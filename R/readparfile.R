readparfile <-
function(parfile=NA, starttime=NA, quiet=FALSE, env=NULL){
  #Procedure to setup Parameter Space (Read .par File) {{{
  if (is.null(env)) { env<-environment() }
  #}}}

  #Initialise Warning Variable {{{
  parwarning<-NULL
  #}}}
  #Check Calling Syntax {{{
  if (is.na(parfile)) {
    stop("Parameter file not supplied. To create the default parameter file, run MeasureFluxes('--makepar').")
  }
  if (is.na(starttime)) {
    parwarning<-c(parwarning,"Start time not supplied - Using current clock time")
    starttime=proc.time()[3]
  }
  #}}}

  #Print Header {{{
  if (!quiet) { cat(paste('------------------------------------------------------\n'))
                cat(paste('   LAMBDAR     version ',packageVersion("LAMBDAR"),': MeasureFluxes\n'))
                cat(paste('------------------------------------------------------\n'))
                cat('Initialising Workspace {\n')
                cat('   Reading Parameter File ') }
  #}}}

  #Test Reading of Parameter File {{{
  no.params<-try(max(count.fields(parfile)),silent=TRUE)
  if (class(no.params)=="try-error") {
  #Stop on Error
    if (grepl("cannot open the connection",no.params[1])) {
      cause="File not found"
    } else {
      cause="Cause unknown"
    }
    stop(paste("Parameter file read failed:",cause))
  } else if (!is.finite(no.params)) {
    cause="Parameter File Empty"
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
    stop("RootDirectory Parameter not in Parameter File")
  }
  #Ensure path ends in a '/'
  if (lastnchar(pathroot,1) != '/') { pathroot<-paste(pathroot,'/',sep="") }
  #}}}

  #Root Directory path {{{
  ID="WorkingDirectory"
  ind<-which(params[ID,]!="")
  pathwork<-params[ID,ind]
  if ((length(ind)==0)||(is.na(pathroot))) {
    stop("WorkingDirectory Parameter not in Parameter File")
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
    stop("OutputDirectory not in Parameter File")
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
    parwarning<-c(parwarning,"BeamArea_SqAS Parameter not present in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"PSFConvolve Parameter table read failed; Using 0 (FALSE)")
        psffilt<-0
      }
      if (is.na(psffilt)) {
        parwarning<-c(parwarning,"PSFConvolve Parameter not in Parameter File; Using 0 (FALSE)")
        psffilt<-0
      }
    } else {
      parwarning<-c(parwarning,"PSFConvolve Parameter not in Parameter File; Using 0 (FALSE)")
      psffilt<-0
    }
  }
  psffilt<-psffilt==1
  #}}}

  #Do we want PSF Matched Photometry? {{{
  ID="PSFWeighted"
  ind<-which(params[ID,]!="")
  PSFWeighted<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(PSFWeighted))) {
    if ((length(ind)==1)) {
      PSFWeighted<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(PSFWeighted)=="try-error") {
        parwarning<-c(parwarning,"PSFWeighted Parameter table read failed; Using 0 (FALSE)")
        PSFWeighted<-0
      }
      if (is.na(PSFWeighted)) {
        parwarning<-c(parwarning,"PSFWeighted Parameter not in Parameter File; Using 0 (FALSE)")
        PSFWeighted<-0
      }
    } else {
      parwarning<-c(parwarning,"PSFWeighted Parameter not in Parameter File; Using 0 (FALSE)")
      PSFWeighted<-0
    }
  }
  PSFWeighted<-(PSFWeighted==1)
  #}}}

  #PSF map filename {{{
  ID="PSFMap"
  ind<-which(params[ID,]!="")
  psfmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(psfmap))) {
    parwarning<-c(parwarning,"PSFMap Parameter not in Paramter File; Using NONE")
    psfmap<-"NONE"
  }
  #Determine if provided psfmap is an image or filelist {{{
  if ((length(psfmap)==1)&(psfmap!="NONE")&(!grepl(".fits", psfmap))) {
    #One file provided without .fits extension - must be filelist
    psfmap<-try(c(t(read.table(file.path(pathroot,psfmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(psfmap)=="try-error") {
      parwarning<-c(parwarning,"PSFMap Parameter table read failed; Using NONE")
      psfmap<-"NONE"
    }
    if (is.na(psfmap)) {
      parwarning<-c(parwarning,"PSFMap Parameter not in Paramter File; Using NONE")
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
          parwarning<-c(parwarning,"Gauss_FWHM_AS Parameter table read failed; Using 0")
          gauss_fwhm_as<-0.0
        }
        if (is.na(gauss_fwhm_as)) {
          parwarning<-c(parwarning,"Gauss_FWHM_AS Parameter not in Parameter File; Using 0")
          gauss_fwhm_as<-0.0
        }
      } else {
        parwarning<-c(parwarning,"Gauss_FWHM_AS Parameter not in Parameter File; Using 0")
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
      str<-paste("Loops with bad parameters:",paste(which(psfmap=="NONE" & gauss_fwhm_as==0.0),collapse=", ",sep=""))
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
  #if (any(PSFWeighted & !psffilt)) {
    #message("WARNING: You've asked for PSF Matched Photometry, and yet specified No PSF Convolution... You could be in for some funky results")
    #parwarning<-c(parwarning,"You've asked for PSF Matched Photometry, and yet specified No PSF Convolution... You could be in for some funky results")
  #}
  #}}}

  #Perform Contaminant removal? {{{
  nocontammap<-NULL
  nNNs<-10
  checkContam<-FALSE
  ID="RemoveContam"
  ind<-which(params[ID,]!="")
  filtcontam<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(filtcontam))) {
    if ((length(ind)==1)) {
      filtcontam<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(filtcontam)=="try-error") {
        parwarning<-c(parwarning,"RemoveContam Parameter table read failed; Using 0 (FALSE)")
        filtcontam<-FALSE
      } else { filtcontam<-(filtcontam==1) }
      if (is.na(filtcontam)) {
        parwarning<-c(parwarning,"RemoveContam Parameter not in Parameter File; Using 0 (FALSE)")
        filtcontam<-FALSE
      }
    } else {
      parwarning<-c(parwarning,"RemoveContam Parameter not in Parameter File; Using 0 (FALSE)")
      filtcontam<-FALSE
    }
  } else { filtcontam<-filtcontam==1 }
  #}}}

  if ( filtcontam ) {
    #Check for irrelevant contaminants? {{{
    ID="CheckContam"
    ind<-which(params[ID,]!="")
    checkContam<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(checkContam))) {
      if ((length(ind)==1)) {
        checkContam<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(checkContam)=="try-error") {
          parwarning<-c(parwarning,"CheckContam Parameter table read failed; Using 0 (FALSE)")
          checkContam<-FALSE
        } else { checkContam<-(checkContam==1) }
        if (is.na(checkContam)) {
          parwarning<-c(parwarning,"CheckContam Parameter not in Parameter File; Using 0 (FALSE)")
          checkContam<-FALSE
        }
      } else {
        parwarning<-c(parwarning,"CheckContam Parameter not in Parameter File; Using 0 (FALSE)")
        checkContam<-FALSE
      }
    } else { checkContam<-checkContam==1 }
    #}}}
    #Check for irrelevant contaminants? {{{
    ID="nNearestCheck"
    ind<-which(params[ID,]!="")
    nNNs<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(nNNs))) {
      if ((length(ind)==1)) {
        nNNs<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(nNNs)=="try-error") {
          parwarning<-c(parwarning,"nNearestCheck Parameter table read failed; Using 10")
          nNNs<-10
        } else { nNNs<-(nNNs==1) }
        if (is.na(nNNs)) {
          parwarning<-c(parwarning,"nNearestCheck Parameter not in Parameter File; Using 10")
          nNNs<-10
        }
      } else {
        parwarning<-c(parwarning,"nNearestCheck Parameter not in Parameter File; Using 10")
        nNNs<-10
      }
    }
    #}}}
    #Contaminant Image Filename {{{
    ID="NoContamImageFile"
    ind<-which(params[ID,]!="")
    nocontammap<-params[ID,ind]
    if ((length(ind)==0)||is.na(nocontammap)) {
      parwarning<-c(parwarning,"NoContamImageFile Parameter not in Parameter File; Using 'NoContamResidualImage.fits'")
      nocontammap<-"NoContamResidualImage.fits"
    }
    #}}}
  }

  #Name of Source Catalogue {{{
  ID="Catalogue"
  ind<-which(params[ID,]!="")
  catalogue<-params[ID,ind]
  if ((length(ind)==0)||(is.na(catalogue))) {
    stop("Catalogue Parameter not in Parameter File")
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
    parwarning<-c(parwarning,"CatIDColumnLabel Parameter not in Parameter File; Using 'CATAID'")
    catalab<-"CATAID"
  }
  #}}}

  #What is the title of the Catalogue's RA Column? {{{
  ID="RAColumnLabel"
  ralab<-params[ID,1]
  if (is.na(ralab)) {
    parwarning<-c(parwarning,"RAColumnLabel Parameter not in Parameter File; Using 'ALPHA_J2000'")
    ralab<-"ALPHA_J2000"
  }#}}}

  #What is the title of the Catalogue's Dec Column? {{{
  ID="DecColumnLabel"
  declab<-params[ID,1]
  if (is.na(declab)) {
    parwarning<-c(parwarning,"DecColumnLabel Parameter not in Parameter File; Using 'DELTA_J2000'")
    declab<-"DELTA_J2000"
  }#}}}

  #What is the title of the Catalogue's Theta Column? {{{
  ID="ThetaColumnLabel"
  thetalab<-params[ID,1]
  if (is.na(thetalab)) {
    parwarning<-c(parwarning,"ThetaColumnLabel Parameter not in Parameter File; Using 'THETA_J2000'")
    thetalab<-"THETA_J2000"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMajColumnLabel"
  semimajlab<-params[ID,1]
  if (is.na(semimajlab)) {
    parwarning<-c(parwarning,"SemiMajColumnLabel Parameter not in Parameter File; Using 'SEMIMAJ_AS'")
    semimajlab<-"SEMIMAJ_AS"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMinColumnLabel"
  semiminlab<-params[ID,1]
  if (is.na(semiminlab)) {
    parwarning<-c(parwarning,"SemiMinColumnLabel Parameter not in Parameter File; Using 'SEMIMIN_AS'")
    semiminlab<-"SEMIMIN_AS"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="ContamColumnLabel"
  contamlab<-params[ID,1]
  if (is.na(contamlab)) {
    parwarning<-c(parwarning,"ContamColumnLabel Parameter not in Parameter File; Using 'CONTAM'")
    contamlab<-"CONTAM"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="FluxWgtColumnLabel"
  fluxweightlab<-params[ID,1]
  if (is.na(fluxweightlab)) {
    parwarning<-c(parwarning,"FluxWgtColumnLabel Parameter not in Parameter File; Using 'FLUXWEIGHT'")
    fluxweightlab<-"FLUXWEIGHT"
  }#}}}

  #Name of Data Image {{{
  ID="DataMap"
  ind<-which(params[ID,]!="")
  datamap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(datamap))) {
    stop("DataMap Parameter not in Parameter File")
  }
  #Determine if provided datamap is an image or filelist {{{
  if ((length(datamap)==1)&(datamap!="NONE")&(!grepl(".fits", datamap))) {
    #One file provided without .fits extension - must be filelist
    datamap<-try(c(t(read.table(file.path(pathroot,datamap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(datamap)=="try-error") {
      #Stop on Error
      stop("Datamap Parameter table read failed")
    }
    if (is.na(datamap)) {
      stop("DataMap Parameter not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Name of Error Map {{{
  ID="ErrorMap"
  ind<-which(params[ID,]!="")
  errormap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(errormap))) {
    parwarning<-c(parwarning,"ErrorMap Parameter not in Parameter File; Using 'NONE'")
    errormap<-"NONE"
  }
  #Determine if provided errormap is an image or filelist {{{
  if ((length(errormap)==1)&(is.na(as.numeric(errormap)))&(errormap!="NONE")&(!grepl(".fits", errormap))) {
    #One file provided without .fits extension - must be filelist
    errormap<-try(c(t(read.table(file.path(pathroot,errormap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(errormap)=="try-error") {
      parwarning<-c(parwarning,"ErrorMap Parameter table read failed; Using 'NONE'")
      errormap<-"NONE"
    }
    if (is.na(errormap)) {
      parwarning<-c(parwarning,"ErrorMap Parameter not in Parameter File; Using 'NONE'")
      errormap<-"NONE"
    }
  }
  #}}}
  #}}}

  #Name of Mask Map {{{
  ID="MaskMap"
  ind<-which(params[ID,]!="")
  maskmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(maskmap))) {
    parwarning<-c(parwarning,"MaskMap Parameter not in Parameter File; Using 'NONE'")
    maskmap<-"NONE"
  }
  #Determine if provided maskmap is an image or filelist {{{
  if ((length(maskmap)==1)&(maskmap!="NONE")&(!grepl(".fits", maskmap))) {
    #One file provided without .fits extension - must be filelist
    maskmap<-try(c(t(read.table(file.path(pathroot,maskmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(maskmap)=="try-error") {
      parwarning<-c(parwarning,"MaskMap Parameter table read failed; Using 'NONE'")
      maskmap<-"NONE"
    }
    if (is.na(maskmap)) {
      parwarning<-c(parwarning,"MaskMap Parameter not in Parameter File; Using 'NONE'")
      maskmap<-"NONE"
    }
  }
  #}}}
  #}}}

  #Name of Weight Map {{{
  ID="WeightMap"
  ind<-which(params[ID,]!="")
  wgtmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(wgtmap))) {
    parwarning<-c(parwarning,"WeightMap Parameter not in Parameter File; Using 'NONE'")
    wgtmap<-"NONE"
  }
  #Determine if provided weightmap is an image or filelist {{{
  if ((length(wgtmap)==1)&(wgtmap!="NONE")&(!grepl(".fits", wgtmap))) {
    #One file provided without .fits extension - must be filelist
    wgtmap<-try(c(t(read.table(file.path(pathroot,wgtmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(wgtmap)=="try-error") {
      parwarning<-c(parwarning,"WeightMap Parameter table read failed; Using 'NONE'")
      wgtmap<-"NONE"
    }
    if (is.na(wgtmap)) {
      parwarning<-c(parwarning,"WeightMap Parameter not in Parameter File; Using 'NONE'")
      wgtmap<-"NONE"
    }
  }
  #}}}
  #}}}

  #Zero Point of Weight Map {{{
  ID="WeightMapZP"
  ind<-which(params[ID,]!="")
  wgtzp<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(wgtzp))) {
    parwarning<-c(parwarning,"WeightMapZP Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"DataExtn Parameter table read failed; Using 0")
        extn<-0
      }
      if (is.na(extn)) {
        parwarning<-c(parwarning,"DataExtn Parameter not in Parameter File; Using 0")
        extn<-0
      }
    } else {
      parwarning<-c(parwarning,"DataExtn Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"ErrorExtn Parameter table read failed; Using 0")
        extnerr<-0
      }
      if (is.na(extnerr)) {
        parwarning<-c(parwarning,"ErrorExtn Parameter not in Parameter File; Using 0")
        extnerr<-0
      }
    } else {
      parwarning<-c(parwarning,"MaskExtn Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"MaskExtn Parameter table read failed; Using 0")
        extnmask<-0
      }
      if (is.na(extnmask)) {
        parwarning<-c(parwarning,"MaskExtn Parameter not in Parameter File; Using 0")
        extnmask<-0
      }
    } else {
      parwarning<-c(parwarning,"MaskExtn Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"WeightExtn Parameter table read failed; Using 0")
        extnwgt<-0
      }
      if (is.na(extnwgt)) {
        parwarning<-c(parwarning,"WeightExtn Parameter not in Parameter File; Using 0")
        extnwgt<-0
      }
    } else {
      parwarning<-c(parwarning,"WeightExtn Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"PointSources Parameter table read failed; Using 0 (FALSE)")
        forcepointsources<-0
      }
      if (is.na(forcepointsources)) {
        parwarning<-c(parwarning,"PointSources Parameter not in Parameter File; Using 0 (FALSE)")
        forcepointsources<-0
      }
    } else {
      parwarning<-c(parwarning,"PointSources Parameter not in Parameter File; Using 0 (FALSE)")
      forcepointsources<-0
    }

  }
  forcepointsources<-(forcepointsources==1)
  #}}}

  #Check for wayward inputs {{{
  #if (any(PSFWeighted & !forcepointsources)) {
    #message("WARNING: You've asked for PSF Matched Photometry, and not forced Point Sources to be used. PSF matched is designed for Point Sources only, so this could behave poorly.")
    #parwarning<-c(parwarning,"You've asked for PSF Matched Photometry, and not forced Point Sources to be used. PSF matched is designed for Point Sources only, so this could behave poorly.")
  #}
  #}}}


  #Error Map scale factor #{{{
  ID="EFactor"
  ind<-which(params[ID,]!="")
  Efactor<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(Efactor))) {
    if ((length(ind)==1)) {
      Efactor<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(Efactor)=='try-error') {
        parwarning<-c(parwarning,"EFactor Parameter table read failed; Using 0")
        Efactor<-1
      }
      if (is.na(Efactor)) {
        parwarning<-c(parwarning,"EFactor Parameter not in Parameter File; Using 0")
        Efactor<-1
      }
    } else {
      parwarning<-c(parwarning,"EFactor Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"FluxCorr Parameter table read failed; Using 1")
        fluxcorr<-1
      }
      if (is.na(fluxcorr)) {
        parwarning<-c(parwarning,"FluxCorr Parameter not in Parameter File; Using 1")
        fluxcorr<-1
      }
    } else {
      parwarning<-c(parwarning,"FluxCorr Parameter not in Parameter File; Using 1")
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
        parwarning<-c(parwarning,"CropImage Parameter table read failed; Using 0 (FALSE)")
        cropimage<-0
      }
      if (is.na(cropimage)) {
        parwarning<-c(parwarning,"CropImage Parameter not in Parameter File; Using 0 (FALSE)")
        cropimage<-0
      }
    } else {
      parwarning<-c(parwarning,"CropImage Parameter not in Parameter File; Using 0 (FALSE)")
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
      parwarning<-c(parwarning,"CropFitsName Parameter not in Parameter File; Using 'croppedimage'")
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
      if (length(ind)==1) {
        #Try Reading table:
        ra0<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(ra0)=="try-error") {
          #Warn on Error
          parwarning<-c(parwarning,"CropImRA0 Parameter table read failed; Using -999")
          ra0<- -999
        }
        if (is.na(ra0)) {
          #Warn on Error
          parwarning<-c(parwarning,"CropImRA0 Parameter not in parameter file; Using -999")
          ra0<- -999
        }
      } else {
        parwarning<-c(parwarning,"CropImRA0 Parameter not in Parameter File; Using -999")
        ra0<- -999
      }
    }
    #}}}
    #Cropped image Dec {{{
    ID="CropImDec0"
    ind<-which(params[ID,]!="")
    dec0<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(dec0))) {
      if (length(ind)==1) {
        #Try Reading table:
        dec0<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(dec0)=="try-error") {
          #Warn on Error
          parwarning<-c(parwarning,"CropImDec0 Parameter table read failed; Using -999")
          dec0<- -999
        }
        if (is.na(dec0)) {
          #Warn on Error
          parwarning<-c(parwarning,"CropImDec0 Parameter not in parameter file; Using -999")
          dec0<- -999
        }
      } else {
        parwarning<-c(parwarning,"CropImDec0 Parameter not in Parameter File; Using -999")
        dec0<- -999
      }
    }
    #}}}
    #Cropped image radius {{{
    ID="CropImRad"
    ind<-which(params[ID,]!="")
    cutrad<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(cutrad))) {
      if (length(ind)==1) {
        #Try Reading table:
        cutrad<-try(as.numeric(c(t(read.table(file.path(pathroot,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(cutrad)=="try-error") {
          #Warn on Error
          parwarning<-c(parwarning,"CropImRad Parameter table read failed; Using 0.5")
          cutrad<- 0.5
        }
        if (is.na(cutrad)) {
          #Warn on Error
          parwarning<-c(parwarning,"CropImRad Parameter not in parameter file; Using 0.5")
          cutrad<- 0.5
        }
      } else {
        parwarning<-c(parwarning,"CropImRad Parameter not in Parameter File; Using 0.5")
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
        parwarning<-c(parwarning,"Confusion_Jy Parameter table read failed; Using 1")
        conf<-1
      }
      if (is.na(conf)) {
        parwarning<-c(parwarning,"Confusion_Jy Parameter not in Parameter File; Using 1")
        conf<-1
      }
    } else {
      parwarning<-c(parwarning,"Confusion_Jy Parameter not in Parameter File; Using 1")
      conf<-1
    }
  }
  #}}}

  #Number of Processors available for computations {{{
  ID="nProcessors"
  ncores<-as.numeric(params[ID,1])
  if (is.na(ncores)) {
    parwarning<-c(parwarning,"nProcessors Parameter not in Parameter File; Using 1")
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
        parwarning<-c(parwarning,"AngularOffset Parameter table read failed; Using 0 (N0E90)")
        angoffset<-0
      }
      if (is.na(angoffset)) {
        parwarning<-c(parwarning,"AngularOffset Parameter not in Parameter File; Using 0 (N0E90)")
        angoffset<-0
      }
    } else {
      parwarning<-c(parwarning,"AngularOffset Parameter not in Parameter File; Using 0 (N0E90)")
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
        parwarning<-c(parwarning,"MapJyPerBeam Parameter table read failed; Using 0 (FALSE)")
        Jybm<-0
      }
      if (is.na(Jybm)) {
        parwarning<-c(parwarning,"MapJyPerBeam Parameter not in Parameter File; Using 0 (FALSE)")
        Jybm<-0
      }
    } else {
      parwarning<-c(parwarning,"MapJyPerBeam Parameter not in Parameter File; Using 0 (FALSE)")
      Jybm<-0
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
        parwarning<-c(parwarning,"SmoothAper Parameter table read failed; Using 1 (TRUE)")
        resampleaperture<-1
      }
      if (is.na(resampleaperture)) {
        parwarning<-c(parwarning,"SmoothAper Parameter not in Parameter File; Using 1 (TRUE)")
        resampleaperture<-1
      }
    } else {
      parwarning<-c(parwarning,"SmoothAper Parameter not in Parameter File; Using 1 (TRUE)")
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
          parwarning<-c(parwarning,"ResamplingRes Parameter table read failed; Using 3")
          upres<-3
        }
        if (is.na(upres)) {
          parwarning<-c(parwarning,"ResamplingRes Parameter not in Parameter File; Using 3")
          upres<-3
        }
      } else {
        parwarning<-c(parwarning,"ResamplingRes Parameter not in Parameter File; Using 3")
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
          parwarning<-c(parwarning,"ResamplingIters Parameter table read failed; Using 5")
          itersteps<-5
        }
        if (is.na(itersteps)) {
          parwarning<-c(parwarning,"ResamplingIters Parameter not in Parameter File; Using 5")
          itersteps<-5
        }
      } else {
        parwarning<-c(parwarning,"ResamplingIters Parameter not in Parameter File; Using 5")
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
        parwarning<-c(parwarning,"PSFConfidence Parameter table read failed; Using 1")
        confidence<-1
      }
      if (is.na(confidence)) {
        parwarning<-c(parwarning,"PSFConfidence Parameter not in Parameter File; Using 1")
        confidence<-1
      }
    } else {
      parwarning<-c(parwarning,"PSFConfidence Parameter not in Parameter File; Using 1")
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
        parwarning<-c(parwarning,"ApStampWidth Parameter table read failed; Using 1.05")
        defbuff<-1.05
      }
      if (is.na(defbuff)) {
        parwarning<-c(parwarning,"ApStampWidth Parameter not in Parameter File; Using 1.05")
        defbuff<-1.05
      }
    } else {
      parwarning<-c(parwarning,"ApStampWidth Parameter not in Parameter File; Using 1.05")
      defbuff<-1.05
    }
  }
  if (defbuff<1) {
    parwarning<-c(parwarning,"ApStampWidth Value is less than or equal to Unity. Value must be strictly > 1. Setting to 1.05")
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
        parwarning<-c(parwarning,"SourceMaskOnly Parameter table read failed; Using 0 (FALSE)")
        sourcemaskonly<-0
      }
      if (is.na(sourcemaskonly)) {
        parwarning<-c(parwarning,"SourceMaskOnly Parameter not in Parameter File; Using 0 (FALSE)")
        sourcemaskonly<-0
      }
    } else {
      parwarning<-c(parwarning,"SourceMaskOnly Parameter not in Parameter File; Using 0 (FALSE)")
      sourcemaskonly<-0
    }
  }
  sourcemaskonly<-(sourcemaskonly==1)
  #}}}

  #Do we want to output the source mask at all?
  if (any(!sourcemaskonly)) {
    ID="WriteSourceMask"
    ind<-which(params[ID,]!="")
    sourcemaskout<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(sourcemaskout)) {
      if (length(ind)==1) {
        sourcemaskout<-try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
        if (class(sourcemaskout)=="try-error") {
          #Warn on Error
          parwarning<-c(parwarning,"WriteSourceMask Parameter table read failed; Using 0 (FALSE)")
          sourcemaskout<-0
        }
        if (is.na(sourcemaskout)) {
          parwarning<-c(parwarning,"WriteSourceMask Parameter not in Parameter File; Using 0 (FALSE)")
          sourcemaskout<-0
        }
      } else {
        parwarning<-c(parwarning,"WriteSourceMask Parameter not in Parameter File; Using 0 (FALSE)")
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
    parwarning<-c(parwarning,"WriteAAMask Parameter not in Parameter File; Using 0 (FALSE)")
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
      parwarning<-c(parwarning,"AllApersFile Parameter not in Parameter File; Using 'AllApertures_Mask.fits'")
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
    parwarning<-c(parwarning,"WriteFAMask Parameter not in Parameter File; Using 0 (FALSE)")
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
      parwarning<-c(parwarning,"ConvApersFile Parameter not in Parameter File; Using 'AllConvolvedApertures_Mask.fits'")
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
    parwarning<-c(parwarning,"WriteDFAMask Parameter not in Parameter File; Using 0 (FALSE)")
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
      parwarning<-c(parwarning,"DeblConvApersFile Parameter not in Parameter File; Using 'AllDeblConvolvedApertures_Mask.fits'")
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
    parwarning<-c(parwarning,"WriteResidMap Parameter not in Parameter File; Using 0 (FALSE)")
    makeresidmap<-FALSE
  } else { makeresidmap<-(makeresidmap==1) }
  #}}}

  #Residual Image filename {{{
  if ( makeresidmap ) {
    #Name of the output Residual image
    ID="ResidImageFile"
    ind<-which(params[ID,]!="")
    residmap<-params[ID,ind]
    if ((length(ind)==0)||is.na(residmap)) {
      parwarning<-c(parwarning,"ResidImageFile Parameter not in Parameter File; Using 'ResidualImage.fits'")
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
    parwarning<-c(parwarning,"WriteTable Parameter not in Parameter File; Using 1 (TRUE)")
    writetab<-TRUE
  } else { writetab<-(writetab==1) }
  #}}}

  #Table Filename {{{
  if ( writetab ) {
    #Name of output Flux table
    ID="TableName"
    ind<-which(params[ID,]!="")
    tableoutname<-params[ID,ind]
    if ((length(ind)==0)||(is.na(tableoutname))) {
      #Warn on Error
      parwarning<-c(parwarning,"TableName Parameter not in Parameter File; Using 'LAMBDAR_Results'")
      tableoutname<-"LAMBDAR_Results"
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
      parwarning<-c(parwarning,"ShowTime Parameter not in Parameter File; Using 0 (FALSE)")
      showtime<-FALSE
    } else { showtime<-(showtime==1) }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Interactive"
  ind<-which(params[ID,]!="")
  interact<-params[ID,ind]
  if ((length(ind)==0)||is.na(interact)) {
    parwarning<-c(parwarning,"Interactive Parameter not in Parameter File; Using 0 (FALSE)")
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
        parwarning<-c(parwarning,"UseMaskLim Parameter table read failed; Using 0.2")
        useMaskLim<-0.2
      }
      if (is.na(useMaskLim)) {
        parwarning<-c(parwarning,"UseMaskLim Parameter not in Parameter File; Using 0.2")
        useMaskLim<-0.2
      }
    } else {
      parwarning<-c(parwarning,"UseMaskLim Parameter not in Parameter File; Using 0.2")
      useMaskLim<-0.2
    }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Diagnostic"
  ind<-which(params[ID,]!="")
  diagnostic<-params[ID,ind]
  if ((length(ind)==0)||is.na(diagnostic)) {
    parwarning<-c(parwarning,"Diagnostic Parameter not in Parameter File; Using 0 (FALSE)")
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
      parwarning<-c(parwarning,"Verbose Parameter not in Parameter File; Using 0 (FALSE)")
      verbose<-FALSE
    } else { verbose<-(verbose==1) }
  }
  #}}}

  #Do we want a sample of the apertures to be output? {{{
  ID="PlotSample"
  ind<-which(params[ID,]!="")
  plotsample<-params[ID,ind]
  if ((length(ind)==0)||is.na(plotsample)) {
    parwarning<-c(parwarning,"PlotSample Parameter not in Parameter File; Using 0 (FALSE)")
    plotsample<-FALSE
  } else { plotsample<-(plotsample==1) }
  ID="PlotAll"
  ind<-which(params[ID,]!="")
  plotall<-params[ID,ind]
  if ((length(ind)==0)||is.na(plotall)) {
    parwarning<-c(parwarning,"PlotAll Parameter not in Parameter File; Using 0 (FALSE)")
    plotall<-FALSE
  } else { plotall<-(plotall==1) }
  if (plotall & !plotsample) {
    parwarning<-c(parwarning,"PlotAll Parameter TRUE but PlotSample Parameter FALSE; Forcing PlotSample Parameter to TRUE")
    plotsample<-TRUE
  }
  #}}}

  #Make Magnitudes in Output? {{{
  ID="Magnitudes"
  ind<-which(params[ID,]!="")
  Magnitudes<-params[ID,ind]
  if ((length(ind)==0)||is.na(Magnitudes)) {
    parwarning<-c(parwarning,"Magnitudes Parameter not present in the Parameter File; Using 1 (TRUE)")
    Magnitudes<-TRUE
  } else { Magnitudes<-(Magnitudes==1) }
  #}}}

  #Magnitude Details {{{
  if (Magnitudes) {
    #AB Vega Magnitude {{{
    ID="ABVegaFlux"
    ind<-which(params[ID,]!="")
    ABvegaflux<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(ABvegaflux))) {
      parwarning<-c(parwarning,"ABVegaFlux Parameter not present in the Parameter File; Using 1.0")
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
          parwarning<-c(parwarning,"MagZeroPoint Parameter table read failed; Using 0.0")
          magZP<-0.0
        }
        if (is.na(magZP)) {
          #Warn on Error
          parwarning<-c(parwarning,"MagZeroPoint Parameter not present/bad in the Parameter File; Using 0.0")
          magZP<-0.0
        }
      } else {
        #Warn on Error
        parwarning<-c(parwarning,"Magnitudes Zero Point not present/bad in the Parameter File; Using 0.0")
        magZP<-0.0
      }
    }
    #}}}

    #Magnitudes Zero Point {{{
    ID="MagZPLabel"
    ind<-which(params[ID,]!="")
    magZPlabel<-params[ID,ind]
    if ((length(ind)==0)||(is.na(magZPlabel))) {
      parwarning<-c(parwarning,"MagZPLabel Parameter not present in the Parameter File; Using 'MagZP'")
      magZPlabel<-"MagZP"
    }
    #}}}
  } else {
    magZP<-NA
    magZPlabel<-"MagZP"
    ABvegaflux<-1.0
  }
  #}}}

  #Saturation Label {{{
  ID="SaturationLabel"
  ind<-which(params[ID,]!="")
  saturlabel<-params[ID,ind]
  if ((length(ind)==0)||(is.na(saturlabel))) {
    parwarning<-c(parwarning,"SaturationLabel Parameter not present in the Parameter File; Using 'SATUR'")
    saturlabel<-"SATUR"
  } else if (length(ind)==1&&(file.exists(file.path(pathroot,saturlabel)))) {
    #One file provided
    saturlabel<-try(c(t(read.table(file.path(pathroot,saturlabel), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(saturlabel)=="try-error") {
      parwarning<-c(parwarning,"SaturationLabel Parameter table read failed; Using 'SATUR'")
      saturlabel<-"SATUR"
    }
    if (is.na(saturlabel)) {
      parwarning<-c(parwarning,"SaturationLabel Parameter not in Parameter File; Using 'SATUR'")
      saturlabel<-"SATUR"
    }
  }
  #}}}

  #Saturation Label {{{
  ID="Saturation"
  ind<-which(params[ID,]!="")
  saturation<-as.numeric(params[ID,ind])
  if ((length(ind)==0)) {
    parwarning<-c(parwarning,"Saturation Parameter not present in the Parameter File; Using Inf")
    saturation<-Inf
  } else if (length(ind)==1&&(file.exists(file.path(pathroot,saturation)))) {
    #One file provided
    saturation<-try(as.numeric(c(t(read.table(file.path(pathroot,saturation), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
    if (class(saturation)=="try-error") {
      parwarning<-c(parwarning,"Saturation Parameter table read failed; Using Inf")
      saturation<-Inf
    }
    if (is.na(saturation)) {
      parwarning<-c(parwarning,"Saturation Parameter not in Parameter File; Using Inf")
      saturation<-Inf
    }
  } else if (any(is.na(saturation))) {
    parwarning<-c(parwarning,"Saturation Parameter is NA in Parameter File; Using Inf")
    saturation[which(is.na(saturation))]<-Inf
  }
  #}}}

  #Gain Label {{{
  ID="GainLabel"
  ind<-which(params[ID,]!="")
  gainlabel<-params[ID,ind]
  if ((length(ind)==0)||(is.na(gainlabel))) {
    parwarning<-c(parwarning,"GainLabel Parameter not present in the Parameter File; Using 'GAIN'")
    gainlabel<-"GAIN"
  } else if (length(ind)==1&&(file.exists(file.path(pathroot,gainlabel)))) {
    #One file provided
    gainlabel<-try(c(t(read.table(file.path(pathroot,gainlabel), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(gainlabel)=="try-error") {
      parwarning<-c(parwarning,"SaturationLabel Parameter table read failed; Using 'SATUR'")
      gainlabel<-"GAIN"
    }
    if (is.na(gainlabel)) {
      parwarning<-c(parwarning,"SaturationLabel Parameter not in Parameter File; Using 'SATUR'")
      gainlabel<-"GAIN"
    }
  }
  #}}}

  #Blanks Correction {{{
  #Do we want to perform a blanks Correction to the errors of object fluxes?
  ID="BlankCor"
  ind<-which(params[ID,]!="")
  BlankCor<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(BlankCor))) {
    if ((length(ind)==1)) {
      BlankCor<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(BlankCor)=="try-error") {
        parwarning<-c(parwarning,"BlankCor Parameter table read failed; Using 0 (FALSE)")
        BlankCor<-0
      }
      if (is.na(BlankCor)) {
        parwarning<-c(parwarning,"BlankCor Parameter not in Parameter File; Using 0 (FALSE)")
        BlankCor<-0
      }
    } else {
      parwarning<-c(parwarning,"BlankCor Parameter not in Parameter File; Using 0 (FALSE)")
      BlankCor<-0
    }
  }
  BlankCor<-BlankCor==1
  #Number of Blanks per Object {{{
  ID="nBlanks"
  ind<-which(params[ID,]!="")
  nBlanks<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(nBlanks))) {
    if ((length(ind)==1)) {
      nBlanks<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(nBlanks)=="try-error") {
        parwarning<-c(parwarning,"nBlanks Parameter table read failed; Using 10")
        nBlanks<-10
      }
      if (is.na(nBlanks)) {
        parwarning<-c(parwarning,"nBlanks Parameter not in Parameter File; Using 10")
        nBlanks<-10
      }
    } else {
      parwarning<-c(parwarning,"nBlanks Parameter not in Parameter File; Using 10")
      nBlanks<-10
    }
  }
  #}}}
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
        parwarning<-c(parwarning,"RanCor Parameter table read failed; Using 0 (FALSE)")
        RanCor<-0
      }
      if (is.na(RanCor)) {
        parwarning<-c(parwarning,"RanCor Parameter not in Parameter File; Using 0 (FALSE)")
        RanCor<-0
      }
    } else {
      parwarning<-c(parwarning,"RanCor Parameter not in Parameter File; Using 0 (FALSE)")
      RanCor<-0
    }
  }
  if(RanCor!="execute") { RanCor<-RanCor==1 }
  #Number of Randoms per Object {{{
  ID="nRandoms"
  ind<-which(params[ID,]!="")
  nRandoms<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(nRandoms))) {
    if ((length(ind)==1)) {
      nRandoms<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(nRandoms)=="try-error") {
        parwarning<-c(parwarning,"nRandoms Parameter table read failed; Using 10")
        nRandoms<-10
      }
      if (is.na(nRandoms)) {
        parwarning<-c(parwarning,"nRandoms Parameter not in Parameter File; Using 10")
        nRandoms<-10
      }
    } else {
      parwarning<-c(parwarning,"nRandoms Parameter not in Parameter File; Using 10")
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
        parwarning<-c(parwarning,"DoSkyEst Parameter table read failed; Using 0 (FALSE)")
        doskyest<-0
      }
      if (is.na(doskyest)) {
        #Warn on Error
        parwarning<-c(parwarning,"DoSkyEst Parameter not in Parameter File; Using 0 (FALSE)")
        doskyest<-0
      }
    } else {
      #Warn on Error
      parwarning<-c(parwarning,"DoSkyEst Parameter not in Parameter File; Using 0 (FALSE)")
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
        parwarning<-c(parwarning,"GetSkyRMS Parameter not in Parameter File; Using 0 (FALSE)")
        getskyrms<-0
      }
      if (is.na(getskyrms)) {
        #Warn on Error
        parwarning<-c(parwarning,"GetSkyRMS Parameter not in Parameter File; Using 0 (FALSE)")
        getskyrms<-0
      }
    } else {
        #Warn on Error
        parwarning<-c(parwarning,"GetSkyRMS Parameter not in Parameter File; Using 0 (FALSE)")
        getskyrms<-0
    }
  }
  getskyrms<-getskyrms==1
  #}}}

  #Sky Estimate Paramters {{{
  if (any(doskyest|getskyrms|BlankCor)) {
    #Sourcemask needed for SkyEstimate. If not TRUE, set to TRUE {{{
    if (any(!sourcemask & (doskyest|getskyrms|BlankCor))) {
      parwarning<-c(parwarning,"Source Mask creation being forced for all runs with doSkyEst/getSkyRMS/BlankCor Parameters TRUE")
      sourcemask<-BlankCor|doskyest|getskyrms|sourcemask
    }
    #}}}

    #Number of iterations used in sky estimation {{{
    ID="SkyEstIters"
    ind<-which(params[ID,]!="")
    skycutiters<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||is.na(skycutiters)) {
      parwarning<-c(parwarning,"SkyEstIters Parameter not in Parameter File; Using 5")
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
          parwarning<-c(parwarning,"SkyEstProbCut Parameter table read failed; Using 3")
          skyprobcut<-3
        }
        if (is.na(skyprobcut)) {
          parwarning<-c(parwarning,"SkyEstProbCut Parameter not in Parameter File; Using 3")
          skyprobcut<-3
        }
      } else {
        parwarning<-c(parwarning,"SkyEstProbCut Parameter not in Parameter File; Using 3")
        skyprobcut<-3
      }
    }
    #}}}

    #Default Sky Value if estimation fails {{{
    ID="SkyDefault"
    ind<-which(params[ID,]!="")
    skydefault<-params[ID,ind]
    if ((length(ind)==0)||is.na(skydefault)) {
      parwarning<-c(parwarning,"SkyDefault Parameter not in Parameter File; Using 'median'")
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
          parwarning<-c(parwarning,"SkyCorrelNoise Parameter table read failed; Using 1")
          correl.noise<-1
        }
        if (is.na(correl.noise)) {
          parwarning<-c(parwarning,"SkyCorrelNoise Parameter not in Parameter File; Using 1")
          correl.noise<-1
        }
      } else {
        parwarning<-c(parwarning,"SkyCorrelNoise Parameter not in Parameter File; Using 1")
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

  #Sourcemask parameters; filename & confidence limit {{{
  smfilename<-NULL
  smConfidenceLim<-NULL
  if ( sourcemask ) {
    if (sourcemaskout) {
      #Name of SourceMask that is output
      ID="SourceMaskFile"
      ind<-which(params[ID,]!="")
      smfilename<-params[ID,ind]
      if (length(ind)==0||is.na(smfilename)||((length(smfilename)==1)&(!grepl(".fits", smfilename)))) {
        if (length(ind)==1) {
          smfilename<-(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
          if (class(smfilename)=="try-error") {
            #Warn on Error
            parwarning<-c(parwarning,"SourceMaskFile Parameter table read failed; Using 'SourceMask.fits'")
            smfilename<-"SourceMask.fits"
          }
          if (is.na(smfilename)) {
            #Warn on Error
            parwarning<-c(parwarning,"SourceMaskFile Parameter not in Parameter File; Using 'SourceMask.fits'")
            smfilename<-"SourceMask.fits"
          }
        } else {
          #Warn on Error
        parwarning<-c(parwarning,"SourceMaskFile Parameter not in Parameter File; Using 'SourceMask.fits'")
          smfilename<-"SourceMask.fits"
        }
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
          parwarning<-c(parwarning,"TransmissionMap Parameter table read failed; Using 0 (FALSE)")
          TransmissionMap<-0
        }
        if (is.na(TransmissionMap)) {
          #Warn on Error
          parwarning<-c(parwarning,"Transmission Map Parameter not in Parameter File; Using 0 (FALSE)")
          TransmissionMap<-0
        }
      } else {
        #Warn on Error
        parwarning<-c(parwarning,"Transmission Map Parameter not in Parameter File; Using 0 (FALSE)")
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
          parwarning<-c(parwarning,"SourceMaskConfLim Parameter table read failed; Using 0.95")
          smConfidenceLim<-0.95
        }
        if (is.na(smConfidenceLim)) {
          #Warn on Error
          parwarning<-c(parwarning,"SourceMaskConfLim Parameter not in Parameter File; Using 0.95")
          smConfidenceLim<-0.95
        }
      } else {
        #Warn on Error
        parwarning<-c(parwarning,"SourceMaskConfLim Parameter not in Parameter File; Using 0.95")
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
        parwarning<-c(parwarning,"MinApRad Parameter not in Parameter File; Using 0")
        MinApRad<-0
      }
      if (is.na(MinApRad)) {
        #Warn on Error
        parwarning<-c(parwarning,"MinApRad Parameter not in Parameter File; Using 0")
        MinApRad<-0
      }
    } else {
      #Warn on Error
      parwarning<-c(parwarning,"MinApRad Parameter not in Parameter File; Using 0")
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
        parwarning<-c(parwarning,"MemorySafe Parameter table read failed; Using 0 (FALSE)")
        memSafe<-0
      }
      if (is.na(memSafe)) {
        memSafe<-0
        parwarning<-c(parwarning,"MemorySafe Parameter not in Parameter File; Using 0 (FALSE)")
      }
    } else {
      memSafe<-0
      parwarning<-c(parwarning,"MemorySafe Parameter not in Parameter File; Using 0 (FALSE)")
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
        parwarning<-c(parwarning,"ApertureConfLimit Parameter table read failed; Using 0.9")
        apLimit<-0.9
      }
      if (is.na(apLimit)) {
        #Warn on Error
        parwarning<-c(parwarning,"ApertureConfLimit Parameter not in Parameter File; Using 0.9")
        apLimit<-0.9
      }
    } else {
      #Warn on Error
      parwarning<-c(parwarning,"ApertureConfLimit Parameter not in Parameter File; Using 0.9")
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
        parwarning<-c(parwarning,"IterateFluxes Parameter table read failed; Using 0 (FALSE)")
        iterateFluxes<-0
      }
      if (is.na(iterateFluxes)) {
        #Warn on Error
        parwarning<-c(parwarning,"IterateFluxes Parameter not in Parameter File; Using 0 (FALSE)")
        iterateFluxes<-0
      }
    } else {
        #Warn on Error
        parwarning<-c(parwarning,"IterateFluxes Parameter not in Parameter File; Using 0 (FALSE)")
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
        parwarning<-c(parwarning,"nIterations Parameter table read failed; Using 2")
        nIterations<-2
      }
      if (is.na(nIterations)) {
        parwarning<-c(parwarning,"nIterations Parameter not in Parameter File; Using 2")
        nIterations<-2
      }
    } else {
      parwarning<-c(parwarning,"nIterations Parameter not in Parameter File; Using 2")
      nIterations<-2
    }
  }
  #}}}

  #What format are the input fluxweights? {{{
  ID="FluxWgtType"
  ind<-which(params[ID,]!="")
  weightType<-tolower(params[ID,ind])
  if ((length(ind)==0)||is.na(weightType)) {
    parwarning<-c(parwarning,"FluxWgtType Parameter not in Parameter File; Using 'scale'")
    weightType<-"scale"
  } else if (weightType!="flux" & weightType!="mag" & weightType!="scale") {
    stop("Fluxweight Type Parameter has unknown value; valid values are 'flux', 'mag', and 'scale'")
  }
  #}}}

  #Do we want to fluxweight Using Pixel Fluxes? {{{
  ID="UsePixelFluxWgts"
  ind<-which(params[ID,]!="")
  usePixelFluxWeights<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(usePixelFluxWeights)) {
    if (length(ind)==1) {
      usePixelFluxWeights<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(usePixelFluxWeights)=="try-error") {
        #Warn on Error
        parwarning<-c(parwarning,"UsePixelFluxWgts Parameter table read failed; Using 0 (FALSE)")
        usePixelFluxWeights<-0
      }
      if (is.na(usePixelFluxWeights)) {
        #Warn on Error
        parwarning<-c(parwarning,"UsePixelFluxWgts Parameter not in Parameter File; Using 0 (FALSE)")
        usePixelFluxWeights<-0
      }
    } else {
      #Warn on Error
      parwarning<-c(parwarning,"UsePixelFluxWgts Parameter not in Parameter File; Using 0 (FALSE)")
      usePixelFluxWeights<-0
    }
  }
  usePixelFluxWeights<-(usePixelFluxWeights==1)
  if (usePixelFluxWeights) { weightType<-"scale" }
  #}}}

  #Do we just want the Deblend Fraction for each aperture? {{{
  ID="GetDeblFrac"
  ind<-which(params[ID,]!="")
  getDeblFrac<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(getDeblFrac)) {
    if (length(ind)==1) {
      getDeblFrac<-as.numeric(try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(getDeblFrac)=="try-error") {
        #Warn on Error
        parwarning<-c(parwarning,"GetDeblFrac Parameter not in Parameter File; Using 0 (FALSE)")
        getDeblFrac<-0
      }
      if (is.na(getDeblFrac)) {
        #Warn on Error
        parwarning<-c(parwarning,"GetDeblFrac Parameter not in Parameter File; Using 0 (FALSE)")
        getDeblFrac<-0
      }
    } else {
        #Warn on Error
        parwarning<-c(parwarning,"GetDeblFrac Parameter not in Parameter File; Using 0 (FALSE)")
        getDeblFrac<-0
    }
  }
  getDeblFrac<-getDeblFrac==1
  #Don't do unwanted things if using this option {{{
  if (any(getDeblFrac & (getskyrms|doskyest))) {
    parwarning<-c(parwarning,"Stopping Sky Estimate for all loops where GetDeblFrac is TRUE, as it is unnecessary")
    if (length(doskyest)==1) {
      doskyest<-!getDeblFrac
    } else if (length(getDeblFrac)==1) {
      doskyest<-FALSE
    } else {
      doskyest[which(getDeblFrac)]<-FALSE
    }
    if (length(getskyrms)==1) {
      getskyrms<-!getDeblFrac
    } else if (length(getDeblFrac)==1) {
      getskyrms<-FALSE
    } else {
      getskyrms[which(getDeblFrac)]<-FALSE
    }
  }
  #}}}
  #}}}

  #Name of Logfile to be output {{{
  ID="LogFile"
  logfile<-params[ID,1]
  if (is.na(logfile)) {
    parwarning<-c(parwarning,"LogFile Parameter not in Parameter File; Using 'LAMBDAR_Log.txt'")
    logfile<-"LAMBDAR_Log.txt"
  }
  #}}}
  #}}}

  # Print any warnings {{{
  if (!is.null(parwarning)) {
    parwarning<-paste(parwarning,collapse="\n     > ")
    cat("{\n    Warnings in Parameter File read:\n     > ")
    cat(parwarning)
    cat("\n   } ")
  }

  #}}}

  # Assign variables to LAMBDAR workspace {{{
  assign("aafilename"       , aafilename       , envir = env) # A
  assign("ABvegaflux"       , ABvegaflux       , envir = env) #
  assign("angoffset"        , angoffset        , envir = env) #
  assign("apLimit"          , apLimit          , envir = env) #
  assign("beamarea_SOM_as"  , beamarea_SOM_as  , envir = env) # B
  assign("BlankCor"         , BlankCor         , envir = env) # B
  assign("conf"             , conf             , envir = env) # C
  assign("checkContam"      , checkContam      , envir = env) #
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
  assign("gainlabel"        , gainlabel        , envir = env) #
  assign("gauss_fwhm_as"    , gauss_fwhm_as    , envir = env) # G
  assign("getDeblFrac"      , getDeblFrac      , envir = env) #
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
  assign("nNNs"             , nNNs             , envir = env) #
  assign("nocontammap"      , nocontammap      , envir = env) #
  assign("ncores"           , ncores           , envir = env) #
  assign("nBlanks"          , nBlanks          , envir = env) #
  assign("nRandoms"         , nRandoms         , envir = env) #
  assign("pathroot"         , pathroot         , envir = env) # P
  assign("pathwork"         , pathwork         , envir = env) #
  assign("pathout"          , pathout          , envir = env) #
  assign("plotsample"       , plotsample       , envir = env) #
  assign("plotall"          , plotall          , envir = env) #
  assign("psfmap"           , psfmap           , envir = env) #
  assign("PSFWeighted"      , PSFWeighted      , envir = env) #
  assign("psffilt"          , psffilt          , envir = env) #
  assign("resampleaperture" , resampleaperture , envir = env) # QR
  assign("ra0"              , ra0              , envir = env) #
  assign("ralab"            , ralab            , envir = env) #
  assign("RanCor"           , RanCor           , envir = env) #
  assign("residmap"         , residmap         , envir = env) #
  assign("saturation"       , saturation       , envir = env) # S
  assign("saturlabel"       , saturlabel       , envir = env) #
  assign("sourcemask"       , sourcemask       , envir = env) #
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
  return=parwarning
  #}}}

}

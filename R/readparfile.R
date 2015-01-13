readparfile <-
function(parfile=NA, starttime=NA, quiet=FALSE, env=NULL){
  #Procedure to setup Parameter Space (Read .par File) {{{
  if (is.null(env)) { env<-environment() }
  #}}}

  #Check Calling Syntax {{{
  if (is.na(parfile)) {
    stop("Parameter file not supplied. To create the default parameter file, run MeasureGamaFluxes('--makepar').")
  }
  if (is.na(starttime)) {
    warning("Start time not supplied - using current clock time")
    starttime=proc.time()[3]
  }
  #}}}

  #Print Header {{{
  if (!quiet) { cat(paste('------------------------------------------------------\n'))
                cat(paste('   LAMBDAR     version ',packageVersion("LAMBDAR"),': MeasureGamaFluxes\n'))
                cat(paste('------------------------------------------------------\n'))
                cat('Initialising Workspace {\n')
                cat('   Reading Parameter File   ') }
  #}}}

  #Test Reading of Parameter File {{{
  params<-try(read.table(parfile, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#", row.names=1, fill=TRUE), silent=TRUE)
  if (class(params)=="try-error") {
    #Stop on Error
    stop("Parameter file read failed")
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
    pathwork<-try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
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
    pathout<-try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
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
    warning("Beamarea Parameter not present in Parameter File")
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
        psffilt<-0
      }
      if (is.na(psffilt)) {
        psffilt<-0
      }
    } else {
      warning("PSF Convolve Flag not in Parameter File")
      psffilt<-0
    }
  }
  nopsf<-(psffilt==0)
  #}}}

  #PSF map filename {{{
  ID="PSFMap"
  ind<-which(params[ID,]!="")
  psfmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(psfmap))) {
    warning("PSF Map Filename not in Paramter File")
    psfmap<-"NONE"
  }
  #Determine if provided psfmap is an image or filelist {{{
  if ((length(psfmap)==1)&(psfmap!="NONE")&(!grepl(".fits", psfmap))) {
    #One file provided without .fits extension - must be filelist
    psfmap<-try(c(t(read.table(file.path(pathroot,psfmap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(psfmap)=="try-error") {
      warning("PSF Map Filename not in Paramter File")
      psfmap<-"NONE"
    }
    if (is.na(psfmap)) {
      warning("PSF Map Filename not in Paramter File")
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
          gauss_fwhm_as<-0.0
        }
        if (is.na(gauss_fwhm_as)) {
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
      print(paste(ind, gauss_fwhm_as[ind],"->", collapse=" : "))
      stop("Parameter file does not provide either PSF map or Gaussian FWHM for One or more files")
    }
    #}}}
  } else {
    #If there is a PSF Map - set gauss fwhm to zero {{{
    gauss_fwhm_as<-0.0
    #}}}
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
        warning("Remove Contaminants Flag not in Parameter File")
        filtcontam<-FALSE
      } else { filtcontam<-(filtcontam==1) }
      if (is.na(filtcontam)) {
        warning("Remove Contaminants Flag not in Parameter File")
        filtcontam<-FALSE
      }
    } else {
      warning("Remove Contaminants Flag not in Parameter File")
      filtcontam<-FALSE
    }
  } else { filtcontam<-(filtcontam==1) }
  #}}}

  #Contaminant Image Filename {{{
  if ( filtcontam ) {
    #Name of the output Residual image
    ID="NoContamImageFile"
    nocontammap<-params[ID,1]
    if (is.na(nocontammap)) {
      warning("Contaminant Subtracted Residual Map Filename not in Parameter File")
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
    warning("Catalogue CATID Column Label not in Parameter File; using 'CATAID'")
    catalab<-"CATAID"
  }
  #}}}

  #What is the title of the Catalogue's RA Column? {{{
  ID="RAColumnLabel"
  ralab<-params[ID,1]
  if (is.na(ralab)) {
    warning("Catalogue RA Column Label not in Parameter File; using 'ALPHA_J2000'")
    ralab<-"ALPHA_J2000"
  }#}}}

  #What is the title of the Catalogue's Dec Column? {{{
  ID="DecColumnLabel"
  declab<-params[ID,1]
  if (is.na(declab)) {
    warning("Catalogue Dec Column Label not in Parameter File; using 'DELTA_J2000'")
    declab<-"DELTA_J2000"
  }#}}}

  #What is the title of the Catalogue's Theta Column? {{{
  ID="ThetaColumnLabel"
  thetalab<-params[ID,1]
  if (is.na(thetalab)) {
    warning("Catalogue Theta Column Label not in Parameter File; using 'THETA_J2000'")
    thetalab<-"THETA_J2000"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMajColumnLabel"
  semimajlab<-params[ID,1]
  if (is.na(semimajlab)) {
    warning("Catalogue SemiMajor Axis Column Label not in Parameter File; using 'SEMIMAJ_AS'")
    semimajlab<-"SEMIMAJ_AS"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMinColumnLabel"
  semiminlab<-params[ID,1]
  if (is.na(semiminlab)) {
    warning("Catalogue SemiMinor Axis Column Label not in Parameter File; using 'SEMIMIN_AS'")
    semiminlab<-"SEMIMIN_AS"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="ContamColumnLabel"
  contamlab<-params[ID,1]
  if (is.na(contamlab)) {
    warning("Catalogue SemiMajor Axis Column Label not in Parameter File; using 'CONTAM'")
    contamlab<-"CONTAM"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="FluxWgtColumnLabel"
  fluxweightlab<-params[ID,1]
  if (is.na(fluxweightlab)) {
    warning("Catalogue FluxWeight Column Label not in Parameter File; using 'FLUXWEIGHT'")
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
    warning("Error Map Path not in Parameter File")
    errormap<-"NONE"
  }
  #Determine if provided errormap is an image or filelist {{{
  if ((length(errormap)==1)&(is.na(as.numeric(errormap)))&(errormap!="NONE")&(!grepl(".fits", errormap))) {
    #One file provided without .fits extension - must be filelist
    errormap<-try(c(t(read.table(file.path(pathroot,errormap), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(errormap)=="try-error") {
      #Stop on Error
      stop("Errormap Filelist read failed")
    }
    if (is.na(errormap)) {
      stop("Data Map Path not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Name of Mask Map {{{
  ID="MaskMap"
  ind<-which(params[ID,]!="")
  maskmap<-params[ID,ind]
  if ((length(ind)==0)||(is.na(maskmap))) {
    warning("Mask Map Path not in Parameter File")
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
    warning("Weight Map Path not in Parameter File")
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
    warning("Weight Map Zero Point not in Parameter File; Using 0")
    wgtzp<-0
  }#}}}

  #Extension number of Data in FITS Header {{{
  ID="DataExtn"
  ind<-which(params[ID,]!="")
  extn<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extn))) {
    warning("FITS Data Extension Value not in Parameter File")
    extn<-0
  }#}}}

  #Extension number of Error Map in FITS Header {{{
  ID="ErrorExtn"
  ind<-which(params[ID,]!="")
  extnerr<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extnerr))) {
    warning("FITS Error Extension Value not in Parameter File")
    extnerr<-0
  }#}}}

  #Extension number of Mask Map in FITS Header {{{
  ID="MaskExtn"
  ind<-which(params[ID,]!="")
  extnmask<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extnmask))) {
    warning("FITS Mask Extension Flag not in Parameter File")
    extnmask<-0
  }#}}}

  #Extension number of Data in FITS Header {{{
  ID="WeightExtn"
  ind<-which(params[ID,]!="")
  extnwgt<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(extnwgt))) {
    warning("FITS Weight Map Extension Value not in Parameter File")
    extnwgt<-0
  }#}}}

  #Do we want to force use of point sources? {{{
  ID="PointSources"
  forcepointsources<-as.numeric(params[ID,1])
  if (is.na(forcepointsources)) {
    warning("Force Point Source Flag not in Parameter File")
    forcepointsources<-FALSE
  } else {
    forcepointsources<-(forcepointsources==1)
  }#}}}

  #Error Map scale factor #{{{
  ID="EFactor"
  ind<-which(params[ID,]!="")
  Efactor<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(Efactor))) {
    warning("Error Map Scale Factor not in Parameter File")
    Efactor<-1
  } #}}}

  #Flux Correction (Scale) Factor {{{
  ID="FluxCorr"
  ind<-which(params[ID,]!="")
  fluxcorr<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(fluxcorr))) {
  warning("Flux Correction Factor not in Parameter File")
  fluxcorr<-1
  }#}}}

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
  cropimage<-params[ID,ind]
  if ((length(ind)==0)||(is.na(cropimage))) {
    warning("Crop Image Flag not in Parameter File")
    cropimage<-FALSE
  } else { cropimage<-(cropimage==1) }
  #}}}

  #Cropped image parameters {{{
  if (cropimage) {
    #What will the cropped image(s) be named {{{
    ID="CropFitsName"
    imfitsoutname<-params[ID,1]
    if (is.na(imfitsoutname)) {
      warning("Cropped Images Filename not in Parameter File")
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
        warning("Cropped Images RA centre not in Parameter File")
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
        warning("Cropped Images Dec Centre not in Parameter File")
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
        warning("Cropped Images cropping Radius not in Parameter File")
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
    warning("Confusion Noise Factor not in Parameter File")
    conf<-0
  }#}}}

  #Number of Processors available for computations {{{
  ID="nProcessors"
  ncores<-as.numeric(params[ID,1])
  if (is.na(ncores)) {
    warning("Number of Processors Value not in Parameter File")
    ncores<-1
  }
  #}}}

  #Angular Offset {{{
  #Is there an offset between the Input Catalogue angles and
  #N0E90 Angular Coordinates?
  ID="AngularOffset"
  angoffset<-params[ID,1]
  if (is.na(angoffset)) {
    warning("Angular Offset Flag not in Parameter File")
    angoffset<-FALSE
  } else { angoffset<-(angoffset==1) }
  #}}}

  #Is the map in Jy per Beam? {{{
  ID="MapJyPerBeam"
  ind<-which(params[ID,]!="")
  Jybm<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(Jybm))) {
    warning("Jansky Per Beam Flag not in Parameter File")
    Jybm<-0
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
        warning("Aperture Resample Flag not in Parameter File")
        resampleaperture<-1
      }
      if (is.na(resampleaperture)) {
        warning("Aperture Resample Flag not in Parameter File")
        resampleaperture<-1
      }
    } else {
      warning("Aperture Resample Flag not in Parameter File")
      resampleaperture<-1
    }
  }
  #}}}

  #Resample Parameters {{{
  if (any(resampleaperture)) {
    #What resolution do we want to upscale by? {{{
    ID="ResamplingRes"
    ind<-which(params[ID,]!="")
    upres<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(upres))) {
      warning("Resampling Resolution Factor not in Parameter File")
      upres<-2
    }
    #}}}
    #How many iterations of upscale do we want? {{{
    ID="ResamplingIters"
    ind<-which(params[ID,]!="")
    itersteps<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(itersteps))) {
      warning("Resampling Iterations Value not in Parameter File")
      itersteps<-10
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
  confidence<-as.numeric(params[ID,1])
  if (is.na(confidence)) {
    warning("PSFConfidence Value not in Parameter File")
    confidence<-0.95
  }  else if (confidence<=0) {
    warning("PSFConfidence Value is less than or equal to zero. Value should be strictly 0<confidence<1. Setting to 0.95")
    confidence<-0.95
  }
  #}}}

  #Size of the aperture stamp as a multiple of the aperture major axis {{{
  ID="ApStampWidth"
  defbuff<-as.numeric(params[ID,1])
  if (is.na(defbuff)) {
    warning("ApStampWidth Value not in Parameter File")
    defbuff<-1.05
  }  else if (defbuff<=1) {
    warning("ApStampWidth Value is less than or equal to Unity. Value must be strictly > 1. Setting to 1.05")
    defbuff<-1.05
  }
  #}}}

  #Do we want to output the source mask only? {{{
  ID="SourceMaskOnly"
  sourcemaskonly<-params[ID,1]
  if (is.na(sourcemaskonly)) {
    warning("Output Source Mask Only Flag not in Parameter File")
    sourcemaskonly<-FALSE
  } else { sourcemaskonly<-(sourcemaskonly==1) }
  #}}}

  #Sourcemask filename {{{
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
  #}}}

  #Do we want to output the All Apertures Mask {{{
  aafilename<-NULL
  ID="WriteAAMask"
  makeaamask<-params[ID,1]
  if (is.na(makeaamask)) {
    warning("Make All Apertures Mask Flag not in Parameter File")
    makeaamask<-FALSE
  } else { makeaamask<-(makeaamask==1) }
  #}}}

  #All Apertures mask file name {{{
  if ( makeaamask ) {
    #Name of the output All Apertures file
    ID="AllApersFile"
    aafilename<-params[ID,1]
    if (is.na(aafilename)) {
      warning("All Apertures Mask Output Filename not in Parameter File")
      aafilename<-"AllApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Convolved Apertures mask {{{
  fafilename<-NULL
  ID="WriteFAMask"
  makefamask<-params[ID,1]
  if (is.na(makefamask)) {
    warning("Make Convolved Apertures Mask Flag not in Parameter File")
    makefamask<-FALSE
  } else { makefamask<-(makefamask==1) }
  #}}}

  #Convolved Apertures Filename {{{
  if ( makefamask ) {
    #Name of the output Convolved Apertures Mask
    ID="ConvApersFile"
    fafilename<-params[ID,1]
    if (is.na(fafilename)) {
      warning("Convolved Apertures Mask Output Filename not in Parameter File")
      fafilename<-"AllConvolvedApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Deblended Convolved Apertures mask {{{
  dfafilename<-NULL
  ID="WriteDFAMask"
  makedfamask<-params[ID,1]
  if (is.na(makedfamask)) {
    warning("Make Deblended Convolved Apertures Mask Flag not in Parameter File")
    makedfamask<-FALSE
  } else { makedfamask<-(makedfamask==1) }
  #}}}

  #Deblended Convolved Apertures Mask filename {{{
  if ( makedfamask ) {
    #Name of the output Convolved Apertures Mask
    ID="DeblConvApersFile"
    dfafilename<-params[ID,1]
    if (is.na(dfafilename)) {
      warning("Convolved Apertures Mask Output Filename not in Parameter File")
      dfafilename<-"AllDeblConvolvedApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Residual image? {{{
  residmap<-NULL
  ID="WriteResidMap"
  makeresidmap<-params[ID,1]
  if (is.na(makeresidmap)) {
    warning("Make Residual Map Flag not in Parameter File")
    makeresidmap<-TRUE
  } else { makeresidmap<-(makeresidmap==1) }
  #}}}

  #Residual Image filename {{{
  if ( makeresidmap ) {
    #Name of the output Residual image
    ID="ResidImageFile"
    residmap<-params[ID,1]
    if (is.na(residmap)) {
      warning("Residual Map Filename not in Parameter File")
      residmap<-"ResidualImage.fits"
    }
  }
  #}}}

  #Do we want to output the Flux table? {{{
  tableoutname<-NULL
  ID="WriteTable"
  writetab<-params[ID,1]
  if (is.na(writetab)) {
    warning("Write Table Flag not in Parameter File")
    writetab<-TRUE
  } else { writetab<-(writetab==1) }
  #}}}

  #Table Filename {{{
  if ( writetab ) {
    #Name of output Flux Table
    ID="TableName"
    tableoutname<-params[ID,1]
    if ((length(ind)==0)||(is.na(tableoutname))) {
        #Warn on Error
        warning("Output Table Filename not in Parameter File")
        tableoutname<-"dfaResults"
      } else {
        if (length(ind)==1) {
        tableoutname<-try(c(t(read.table(file.path(pathroot,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
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
    showtime<-params[ID,1]
    if (is.na(showtime)) {
      warning("ShowTime Flag not in Parameter File")
      showtime<-FALSE
    } else { showtime<-(showtime==1) }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Interactive"
  interact<-params[ID,1]
  if (is.na(interact)) {
    warning("Interactive Flag not in Parameter File")
    interact<-FALSE
  } else { interact<-(interact==1) }
  #}}}

  #What limit do we want for the use of masks what cross the mask edges {{{
  ID="UseMaskLim"
  useMaskLim<-as.numeric(params[ID,1])
  if (is.na(useMaskLim)) {
    warning("UseMaskLim Value not in Parameter File")
    useMaskLim<-0.95
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Diagnostic"
  diagnostic<-params[ID,1]
  if (is.na(diagnostic)) {
    warning("Diagnostic Flag not in Parameter File")
    diagnostic<-FALSE
  } else { diagnostic<-(diagnostic==1) }
  #}}}

  #Do we want Verbose Output in Log File {{{
  if (diagnostic) {
    verbose<-TRUE
  } else {
    ID="Verbose"
    verbose<-params[ID,1]
    if (is.na(verbose)) {
      warning("Verbose Flag not in Parameter File")
      verbose<-FALSE
    } else { verbose<-(verbose==1) }
  }
  #}}}

  #Do we want a sample of the apertures to be output? {{{
  ID="PlotSample"
  plotsample<-params[ID,1]
  if (is.na(plotsample)) {
    warning("Plot Sample Flag not in Parameter File")
    plotsample<-FALSE
  } else { plotsample<-(plotsample==1) }
  #}}}

  #Make Magnitudes in Output? {{{
  ID="Magnitudes"
  Magnitudes<-params[ID,1]
  if (is.na(Magnitudes)) {
    warning("Make Magnitudes Flag not present in the Parameter File")
    Magnitudes<-TRUE
  } else { Magnitudes<-(Magnitudes==1) }
  #}}}

  #AB Vega Magnitude {{{
  ID="ABVegaFlux"
  ind<-which(params[ID,]!="")
  ABvegaflux<-as.numeric(params[ID,ind])
  if (is.na(ABvegaflux)) {
    warning("AB Vega Flux Value not present in the Parameter File; using 1.0")
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
        warning("Magnitudes Zero Point not present/bad in the Parameter File; using 0.0")
        magZP<-0.0
      }
      if (is.na(magZP)) {
        #Warn on Error
        warning("Magnitudes Zero Point not present/bad in the Parameter File; using 0.0")
        magZP<-0.0
      }
    } else {
      #Warn on Error
      warning("Magnitudes Zero Point not present/bad in the Parameter File; using 0.0")
      magZP<-0.0
    }
  }
  #}}}

  #Magnitudes Zero Point {{{
  ID="MagZPLabel"
  ind<-which(params[ID,]!="")
  magZPlabel<-params[ID,ind]
  if ((length(ind)==0)||(is.na(magZPlabel))) {
    warning("FITS Zero Point Label not present in the Parameter File")
    magZPlabel<-"MagZP"
  }
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
        warning("Sky Estimate Flag not present in the Parameter File")
        doskyest<-0
      }
      if (is.na(doskyest)) {
        #Warn on Error
        warning("Sky Estimate Flag not present in the Parameter File")
        doskyest<-0
      }
    } else {
        #Warn on Error
        warning("Sky Estimate Flag not present in the Parameter File")
        doskyest<-0
    }
  }
  doskyest<-(doskyest==1)
  #}}}

  #Calculate the Sky RMS? {{{
  ID="GetSkyRMS"
  getskyrms<-params[ID,1]
  if (is.na(getskyrms)) {
    warning("Get Sky RMS Flag not present in the Parameter File")
    getskyrms<-TRUE
  } else { getskyrms<-(getskyrms==1) }
  #}}}

  #Sky Estimate Paramters {{{
  if (any(doskyest|getskyrms)) {
    #Sourcemask needed for SkyEstimate. If not TRUE, set to TRUE {{{
    if (any(!sourcemask)) {
      warning("Source Mask creation being forced for all runs with Sky Estimate Flag TRUE")
      sourcemask<-doskyest|getskyrms
    }
    #}}}

    #Number of iterations used in sky estimation {{{
    ID="SkyEstIters"
    skycutiters<-params[ID,1]
    if (is.na(skycutiters)) {
      warning("Sky Estimate Iterations Value not present in the Parameter File")
      skycutiters<-5
    }
    #}}}

    #Sigma level used in sky cut {{{
    ID="SkyEstProbCut"
    skyprobcut<-as.numeric(params[ID,1])
    if (is.na(skyprobcut)) {
      warning("Sky Estimate Simga Cut Level not present in the Parameter File")
      skyprobcut<-3
    }
    #}}}

    #Default Sky Value if estimation fails {{{
    ID="SkyDefault"
    skydefault<-params[ID,1]
    if (is.na(skydefault)) {
      warning("Sky Default Value not present in the Parameter File; using median")
      skydefault<-"median"
    }
    #}}}

    #Level of Corellated noise in image {{{
    ID="SkyCorrelNoise"
    ind<-which(params[ID,]!="")
    correl.noise<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(correl.noise))) {
      warning("Sky Correlated Noise Level not present in the Parameter File")
      correl.noise<-1
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
        warning("Minimum Aperture Specification not in Parameter File")
        MinApRad<-0
      }
      if (is.na(MinApRad)) {
        #Warn on Error
        warning("Minimum Aperture Specification not in Parameter File")
        MinApRad<-0
      }
    } else {
      #Warn on Error
      warning("Minimum Aperture Specification not in Parameter File")
      MinApRad<-0
    }
  }
  #}}}

  #Do we want a Memory-Safe run? {{{
  ID="MemorySafe"
  memSafe<-as.numeric(params[ID,1])
  if (is.na(memSafe)) {
   warning("Memory Safe Flag not present in the Parameter File")
   memSafe<-TRUE
  } else { memSafe<-(memSafe==1) }
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
        warning("Aperture Confidence Limit Value not in Parameter File")
        apLimit<-0.9
      }
      if (is.na(apLimit)) {
        #Warn on Error
        warning("Aperture Confidence Limit Value not in Parameter File")
        apLimit<-0.9
      }
    } else {
      #Warn on Error
      warning("Aperture Confidence Limit Value not in Parameter File")
      apLimit<-0.9
    }
  }
  #}}}

  #Do we want to iterate the fluxes? {{{
  ID="IterateFluxes"
  iterateFluxes<-params[ID,1]
  if (is.na(iterateFluxes)) {
    warning("Iterate Fluxes Parameter not present in the Parameter File")
    iterateFluxes<-FALSE
  } else { iterateFluxes<-(iterateFluxes==1) }
  #}}}

  #If Iterating, how many iterations do we want? {{{
  ID="nIterations"
  nIterations<-as.numeric(params[ID,1])
  if (is.na(nIterations)) {
    warning("Number of Iterations Parameter not present in the Parameter File")
    nIterations<-2
  }
  #}}}

  #What format are the input fluxweights? {{{
  ID="FluxWgtType"
  weightType<-tolower(params[ID,1])
  if (is.na(weightType)) {
    warning("Fluxweight Type Parameter not present in the Parameter File")
    weightType<-"scale"
  } else if (weightType!="flux" & weightType!="mag" & weightType!="scale") {
    stop("Fluxweight Type Parameter has unknown value; valid values are 'flux', 'mag', and 'scale'")
  }
  #}}}

  #Do we want to fluxweight using Pixel Fluxes? {{{
  ID="UsePixelFluxWgts"
  usePixelFluxWeights<-params[ID,1]
  if (is.na(usePixelFluxWeights)) {
    warning("Use PixelFlux for Fluxweights Parameter not present in the Parameter File")
    usePixelFluxWeights<-TRUE
  } else { usePixelFluxWeights<-(usePixelFluxWeights==1) }
  if (usePixelFluxWeights) { weightType<-"scale" }
  #}}}

  #Do we want to Create a Simulated Image? {{{
  ID="CreateSimulatedImage"
  createSimImage<-params[ID,1]
  if (is.na(createSimImage)) {
    createSimImage<-FALSE
  } else { createSimImage<-(createSimImage==1) }
  #}}}

  #Name of Logfile to be output {{{
  ID="LogFile"
  logfile<-params[ID,1]
  if (is.na(logfile)) {
    warning("Output Log Filename not present in the Parameter File")
    logfile<-"LAMDBAR_Log.txt"
  }
  #}}}
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
  assign("createSimImage"   , createSimImage   , envir = env) #
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
  assign("pathroot"         , pathroot         , envir = env) # P
  assign("pathwork"         , pathwork         , envir = env) #
  assign("pathout"          , pathout          , envir = env) #
  assign("plotsample"       , plotsample       , envir = env) #
  assign("psfmap"           , psfmap           , envir = env) #
  assign("resampleaperture" , resampleaperture , envir = env) # QR
  assign("ra0"              , ra0              , envir = env) #
  assign("ralab"            , ralab            , envir = env) #
  assign("residmap"         , residmap         , envir = env) #
  assign("sourcemask"       , sourcemask       , envir = env) # S
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

read.par.file <-
function(par.file=NA, start.time=NA, quiet=FALSE, env=NULL){
  #Procedure to setup Parameter Space (Read .par File) {{{
  if (is.null(env)) { env<-environment() }
  #}}}

  #Initialise Warning Variable {{{
  param.warnings<-NULL
  #}}}
  #Check Calling Syntax {{{
  if (is.na(par.file)) {
    stop("Parameter file not supplied. To create the default parameter file, run measure.fluxes('--makepar').")
  }
  if (is.na(start.time)) {
    param.warnings<-c(param.warnings,"Start time not supplied - Using current clock time")
    start.time=proc.time()[3]
  }
  #}}}

  #Print Header {{{
  if (!quiet) { calls<-sys.status()$sys.calls
                cat(paste('------------------------------------------------------\n'))
                cat(paste('   LAMBDAR : version ',packageVersion("LAMBDAR"),': ',strsplit(paste(calls[[length(calls)-2]]),"(",fixed=TRUE)[[1]][1],'\n'))
                cat(paste('------------------------------------------------------\n'))
                cat('Initialising Workspace {\n')
                cat('   Reading Parameter File ') }
  #}}}

  #Test Reading of Parameter File {{{
  no.params<-try(max(count.fields(par.file)),silent=TRUE)
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
  params<-try(read.table(par.file, strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#", row.names=1, fill=TRUE, col.names=1:no.params), silent=TRUE)
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
  path.root<-params[ID,1]
  if ((length(path.root)==0)||(is.na(path.root))) {
    stop("RootDirectory Parameter not in Parameter File")
  }
  #Ensure path ends in a '/'
  if (last.n.char(path.root,1) != '/') { path.root<-paste(path.root,'/',sep="") }
  #}}}

  #Root Directory path {{{
  ID="WorkingDirectory"
  ind<-which(params[ID,]!="")
  path.work<-params[ID,ind]
  if ((length(ind)==0)||(is.na(path.root))) {
    stop("WorkingDirectory Parameter not in Parameter File")
  } else {
    path.work<-try(suppressWarnings(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
    if (class(path.work)=="try-error") {
      path.work<-params[ID,ind]
    }
  }
  #Ensure path ends in a '/'
  if (any(last.n.char(path.work,1) != '/')) { path.work<-paste(path.work,'/',sep="") }
  #}}}

  #Output Directory path {{{
  ID="OutputDirectory"
  ind<-which(params[ID,]!="")
  path.out<-params[ID,ind]
  if ((length(ind)==0)||(is.na(path.out))) {
    stop("OutputDirectory not in Parameter File")
  } else {
    path.out<-try(suppressWarnings(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
    if (class(path.out)=="try-error") {
      path.out<-params[ID,ind]
    }
  }
  #Ensure path ends in a '/'
  if (last.n.char(path.out,1) != '/') { path.out<-paste(path.out,'/',sep="") }
  #}}}

  #Beam area in square arcsec {{{
  ID="BeamArea_SqAS"
  ind<-which(params[ID,]!="")
  beam.area.input.as<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(beam.area.input.as))) {
    param.warnings<-c(param.warnings,"BeamArea_SqAS Parameter not present in Parameter File; Using 0")
    beam.area.input.as<-0
  }
  #}}}

  #Do we want to Convolve the apertures with a PSF {{{
  ID="PSFConvolve"
  ind<-which(params[ID,]!="")
  psf.filt<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(psf.filt))) {
    if ((length(ind)==1)) {
      psf.filt<-try(as.numeric(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
      if (class(psf.filt)=="try-error") {
        param.warnings<-c(param.warnings,"PSFConvolve Parameter table read failed; Using 0 (FALSE)")
        psf.filt<-0
      }
      if (is.na(psf.filt)) {
        param.warnings<-c(param.warnings,"PSFConvolve Parameter not in Parameter File; Using 0 (FALSE)")
        psf.filt<-0
      }
    } else {
      param.warnings<-c(param.warnings,"PSFConvolve Parameter not in Parameter File; Using 0 (FALSE)")
      psf.filt<-0
    }
  }
  psf.filt<-psf.filt==1
  #}}}

  #Do we want PSF Matched Photometry? {{{
  ID="PSFWeighted"
  ind<-which(params[ID,]!="")
  psf.weighted<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(psf.weighted))) {
    if ((length(ind)==1)) {
      psf.weighted<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(psf.weighted)=="try-error") {
        param.warnings<-c(param.warnings,"PSFWeighted Parameter table read failed; Using 0 (FALSE)")
        psf.weighted<-0
      }
      if (is.na(psf.weighted)) {
        param.warnings<-c(param.warnings,"PSFWeighted Parameter not in Parameter File; Using 0 (FALSE)")
        psf.weighted<-0
      }
    } else {
      param.warnings<-c(param.warnings,"PSFWeighted Parameter not in Parameter File; Using 0 (FALSE)")
      psf.weighted<-0
    }
  }
  optimal.aper<-psf.weighted==-1
  psf.weighted<-(psf.weighted==1)
  #}}}

  #PSF map filename {{{
  ID="PSFMap"
  ind<-which(params[ID,]!="")
  psf.map<-params[ID,ind]
  if ((length(ind)==0)||(is.na(psf.map))) {
    param.warnings<-c(param.warnings,"PSFMap Parameter not in Paramter File; Using NONE")
    psf.map<-"NONE"
  }
  #Determine if provided psf.map is an image or filelist {{{
  if ((length(psf.map)==1)&(psf.map!="NONE")&(!grepl(".fits", psf.map,ignore.case=TRUE))) {
    #One file provided without .fits extension - must be filelist
    psf.map<-try(c(t(read.table(file.path(path.root,psf.map), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(psf.map)=="try-error") {
      param.warnings<-c(param.warnings,"PSFMap Parameter table read failed; Using NONE")
      psf.map<-"NONE"
    }
    if (is.na(psf.map)) {
      param.warnings<-c(param.warnings,"PSFMap Parameter not in Paramter File; Using NONE")
      psf.map<-"NONE"
    }
  }
  #}}}
  #}}}

  #If no PSF map, get gaussian FWHM {{{
  if (any(psf.map=="NONE")) {
    #FWHM of seeing gaussian
    ID="Gauss_FWHM_AS"
    ind<-which(params[ID,]!="")
    gauss.fwhm.arcsec<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(gauss.fwhm.arcsec))) {
      if ((length(ind)==1)) {
        gauss.fwhm.arcsec<-try(as.numeric(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
        if (class(gauss.fwhm.arcsec)=="try-error") {
          param.warnings<-c(param.warnings,"Gauss_FWHM_AS Parameter table read failed; Using 0")
          gauss.fwhm.arcsec<-0.0
        }
        if (is.na(gauss.fwhm.arcsec)) {
          param.warnings<-c(param.warnings,"Gauss_FWHM_AS Parameter not in Parameter File; Using 0")
          gauss.fwhm.arcsec<-0.0
        }
      } else {
        param.warnings<-c(param.warnings,"Gauss_FWHM_AS Parameter not in Parameter File; Using 0")
        gauss.fwhm.arcsec<-0.0
      }
    }
    #Make sure PSF maps and Gauss FWHM vals are conformable {{{
    if (length(gauss.fwhm.arcsec)!=length(psf.map)) {
      gauss.fwhm.arcsec<-rep(gauss.fwhm.arcsec[1], length(psf.map))
    }#}}}

    #Make sure files with PSF maps have Gauss FWHM vals set to 0 {{{
    ind<-which(psf.map!="NONE")
    if (length(ind)>0) { gauss.fwhm.arcsec[ind]<-0.0 }
    #}}}

    #If we want convolution, there is no PSF, and no gaussian FWHM provided - ERROR {{{
    ind<-which(psf.map=="NONE")
    if ((psf.filt)&(any(gauss.fwhm.arcsec[ind]==0.0))) {
      cat(" - Error\n")
      str<-paste("Loops with bad parameters:",paste(which(psf.map=="NONE" & gauss.fwhm.arcsec==0.0),collapse=", ",sep=""))
      stop(paste("Parameter file does not provide either PSF map or Gaussian FWHM for One or more files.\n",str,sep=""))
    }
    #}}}
  } else {
    #If there is a PSF Map - set gauss fwhm to zero {{{
    gauss.fwhm.arcsec<-0.0
    #}}}
  }
  #}}}

  #PSF FWHM Label {{{
  ID="PSFLabel"
  ind<-which(params[ID,]!="")
  psf.label<-params[ID,ind]
  if ((length(ind)==0)||(is.na(psf.label))) {
    param.warnings<-c(param.warnings,"PSFLabel Parameter not present in the Parameter File; Using 'PSFSEE'")
    psf.label<-"PSFSEE"
  } else if (length(ind)==1&&(file.exists(file.path(path.root,psf.label)))) {
    #One file provided
    psf.label<-try(c(t(read.table(file.path(path.root,psf.label), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(psf.label)=="try-error") {
      param.warnings<-c(param.warnings,"psfLabel Parameter table read failed; Using 'PSFSEE'")
      psf.label<-"PSFSEE"
    }
    if (is.na(psf.label)) {
      param.warnings<-c(param.warnings,"psfLabel Parameter not in Parameter File; Using 'PSFSEE'")
      psf.label<-"PSFSEE"
    }
  }
  #}}}

  #PSF FWHM Label Type {{{
  ID="PSFLabelType"
  ind<-which(params[ID,]!="")
  psf.label.type<-params[ID,ind]
  if ((length(ind)==0)||(is.na(psf.label.type))) {
    param.warnings<-c(param.warnings,"PSFLabelType Parameter not present in the Parameter File; Using 'FWHM'")
    psf.label.type<-"FWHM"
  } else if (length(ind)==1&&(file.exists(file.path(path.root,psf.label.type)))) {
    #One file provided
    psf.label.type<-try(c(t(read.table(file.path(path.root,psf.label.type), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(psf.label.type)=="try-error") {
      param.warnings<-c(param.warnings,"psfLabel Parameter table read failed; Using 'PSFSEE'")
      psf.label.type<-"FWHM"
    }
    if (is.na(psf.label.type)) {
      param.warnings<-c(param.warnings,"psfLabel Parameter not in Parameter File; Using 'PSFSEE'")
      psf.label.type<-"FWHM"
    }
  }
  #}}}

  #Flag loops with no provided PSF {{{
  no.psf<-((psf.map=="NONE") & (gauss.fwhm.arcsec==0.0))
  #}}}

  #Perform Contaminant removal? {{{
  no.contam.map<-NULL
  num.nearest.neighbours<-10
  check.contam<-FALSE
  ID="RemoveContam"
  ind<-which(params[ID,]!="")
  filt.contam<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(filt.contam))) {
    if ((length(ind)==1)) {
      filt.contam<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(filt.contam)=="try-error") {
        param.warnings<-c(param.warnings,"RemoveContam Parameter table read failed; Using 0 (FALSE)")
        filt.contam<-FALSE
      } else { filt.contam<-(filt.contam==1) }
      if (is.na(filt.contam)) {
        param.warnings<-c(param.warnings,"RemoveContam Parameter not in Parameter File; Using 0 (FALSE)")
        filt.contam<-FALSE
      }
    } else {
      param.warnings<-c(param.warnings,"RemoveContam Parameter not in Parameter File; Using 0 (FALSE)")
      filt.contam<-FALSE
    }
  } else { filt.contam<-filt.contam==1 }
  #}}}

  if ( filt.contam ) {
    #Check for irrelevant contaminants? {{{
    ID="CheckContam"
    ind<-which(params[ID,]!="")
    check.contam<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(check.contam))) {
      if ((length(ind)==1)) {
        check.contam<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(check.contam)=="try-error") {
          param.warnings<-c(param.warnings,"CheckContam Parameter table read failed; Using 0 (FALSE)")
          check.contam<-FALSE
        } else { check.contam<-(check.contam==1) }
        if (is.na(check.contam)) {
          param.warnings<-c(param.warnings,"CheckContam Parameter not in Parameter File; Using 0 (FALSE)")
          check.contam<-FALSE
        }
      } else {
        param.warnings<-c(param.warnings,"CheckContam Parameter not in Parameter File; Using 0 (FALSE)")
        check.contam<-FALSE
      }
    } else { check.contam<-check.contam==1 }
    #}}}
    #Number of nearest Neighbours to check? {{{
    ID="nNearestCheck"
    ind<-which(params[ID,]!="")
    num.nearest.neighbours<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(num.nearest.neighbours))) {
      if ((length(ind)==1)) {
        num.nearest.neighbours<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(num.nearest.neighbours)=="try-error") {
          param.warnings<-c(param.warnings,"nNearestCheck Parameter table read failed; Using 10")
          num.nearest.neighbours<-10
        } else { num.nearest.neighbours<-(num.nearest.neighbours==1) }
        if (is.na(num.nearest.neighbours)) {
          param.warnings<-c(param.warnings,"nNearestCheck Parameter not in Parameter File; Using 10")
          num.nearest.neighbours<-10
        }
      } else {
        param.warnings<-c(param.warnings,"nNearestCheck Parameter not in Parameter File; Using 10")
        num.nearest.neighbours<-10
      }
    }
    #}}}
    #Contaminant Image Filename {{{
    ID="NoContamImageFile"
    ind<-which(params[ID,]!="")
    no.contam.map<-params[ID,ind]
    if ((length(ind)==0)||is.na(no.contam.map)) {
      param.warnings<-c(param.warnings,"NoContamImageFile Parameter not in Parameter File; Using 'NoContamResidualImage.fits'")
      no.contam.map<-"NoContamResidualImage.fits"
    }
    #}}}
  }

  #Check for irrelevant contaminants? {{{
  ID="GroupWeights"
  ind<-which(params[ID,]!="")
  group.weights<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(group.weights))) {
    if ((length(ind)==1)) {
      group.weights<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(group.weights)=="try-error") {
        param.warnings<-c(param.warnings,"GroupWeights Parameter table read failed; Using 0 (FALSE)")
        group.weights<-FALSE
      } else { group.weights<-(group.weights==1) }
      if (is.na(group.weights)) {
        param.warnings<-c(param.warnings,"GroupWeights Parameter not in Parameter File; Using 0 (FALSE)")
        group.weights<-FALSE
      }
    } else {
      param.warnings<-c(param.warnings,"GroupWeights Parameter not in Parameter File; Using 0 (FALSE)")
      group.weights<-FALSE
    }
  } else { group.weights<-group.weights==1 }
  #}}}

  #Name of Source Catalogue {{{
  ID="Catalogue"
  ind<-which(params[ID,]!="")
  catalogue<-params[ID,ind]
  if ((length(ind)==0)||(is.na(catalogue))) {
    stop("Catalogue Parameter not in Parameter File")
  }
  #Determine if provided weightmap is an image or filelist {{{
  if ((length(catalogue)==1)&(!grepl(".csv",catalogue,ignore.case=TRUE))&(!grepl(".Rdata",catalogue,ignore.case=TRUE))&(!grepl(".fits", catalogue,ignore.case=TRUE))) {
    #One file provided without .fits extension - must be filelist
    catalogue<-try(c(t(read.table(file.path(path.root,catalogue), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
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
  cata.lab<-params[ID,1]
  if (is.na(cata.lab)) {
    param.warnings<-c(param.warnings,"CatIDColumnLabel Parameter not in Parameter File; Using 'CATAID'")
    cata.lab<-"CATAID"
  }
  #}}}

  #What is the title of the Catalogue's RA Column? {{{
  ID="RAColumnLabel"
  ra.lab<-params[ID,1]
  if (is.na(ra.lab)) {
    param.warnings<-c(param.warnings,"RAColumnLabel Parameter not in Parameter File; Using 'ALPHA_J2000'")
    ra.lab<-"ALPHA_J2000"
  }#}}}

  #What is the title of the Catalogue's Dec Column? {{{
  ID="DecColumnLabel"
  dec.lab<-params[ID,1]
  if (is.na(dec.lab)) {
    param.warnings<-c(param.warnings,"DecColumnLabel Parameter not in Parameter File; Using 'DELTA_J2000'")
    dec.lab<-"DELTA_J2000"
  }#}}}

  #What is the title of the Catalogue's Theta Column? {{{
  ID="ThetaColumnLabel"
  theta.lab<-params[ID,1]
  if (is.na(theta.lab)) {
    param.warnings<-c(param.warnings,"ThetaColumnLabel Parameter not in Parameter File; Using 'THETA_J2000'")
    theta.lab<-"THETA_J2000"
  }#}}}

  #What is the title of the Catalogue's SemiMaj Axis Column? {{{
  ID="SemiMajColumnLabel"
  semimaj.lab<-params[ID,1]
  if (is.na(semimaj.lab)) {
    param.warnings<-c(param.warnings,"SemiMajColumnLabel Parameter not in Parameter File; Using 'SEMIMAJ.arcsec'")
    semimaj.lab<-"SEMIMAJ.arcsec"
  }#}}}

  #What is the title of the Catalogue's SemiMin Axis Column? {{{
  ID="SemiMinColumnLabel"
  semimin.lab<-params[ID,1]
  if (is.na(semimin.lab)) {
    param.warnings<-c(param.warnings,"SemiMinColumnLabel Parameter not in Parameter File; Using 'SEMIMIN.arcsec'")
    semimin.lab<-"SEMIMIN.arcsec"
  }#}}}

  #What is the title of the Catalogue's Contaminant Axis Column? {{{
  ID="ContamColumnLabel"
  contam.lab<-params[ID,1]
  if (is.na(contam.lab)) {
    param.warnings<-c(param.warnings,"ContamColumnLabel Parameter not in Parameter File; Using 'CONTAM'")
    contam.lab<-"CONTAM"
  }#}}}

  #What is the title of the Catalogue's FluxWeight Axis Column? {{{
  ID="FluxWgtColumnLabel"
  flux.weight.lab<-params[ID,1]
  if (is.na(flux.weight.lab)) {
    param.warnings<-c(param.warnings,"FluxWgtColumnLabel Parameter not in Parameter File; Using 'FLUXWEIGHT'")
    flux.weight.lab<-"FLUXWEIGHT"
  }#}}}

  #What is the title of the Catalogue's Grouping Axis Column? {{{
  ID="GroupColumnLabel"
  group.lab<-params[ID,1]
  if (is.na(group.lab)) {
    param.warnings<-c(param.warnings,"GroupColumnLabel Parameter not in Parameter File; Using 'GROUP'")
    group.lab<-"GROUP"
  }#}}}

  #Name of Data Image {{{
  ID="DataMap"
  ind<-which(params[ID,]!="")
  data.map<-params[ID,ind]
  if ((length(ind)==0)||(is.na(data.map))) {
    stop("DataMap Parameter not in Parameter File")
  }
  #Determine if provided data.map is an image or filelist {{{
  if ((length(data.map)==1)&(data.map!="NONE")&(!grepl(".fits", data.map,ignore.case=TRUE))) {
    #One file provided without .fits extension - must be filelist
    data.map<-try(c(t(read.table(file.path(path.root,data.map), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(data.map)=="try-error") {
      #Stop on Error
      stop("Datamap Parameter table read failed")
    }
    if (is.na(data.map)) {
      stop("DataMap Parameter not in Parameter File")
    }
  }
  #}}}
  #}}}

  #Name of Error Map {{{
  ID="ErrorMap"
  ind<-which(params[ID,]!="")
  error.map<-params[ID,ind]
  if ((length(ind)==0)||(is.na(error.map))) {
    param.warnings<-c(param.warnings,"ErrorMap Parameter not in Parameter File; Using 'NONE'")
    error.map<-"NONE"
  }
  #Determine if provided error.map is an image or filelist {{{
  if ((length(error.map)==1)&(is.na(as.numeric(error.map)))&(error.map!="NONE")&(!grepl(".fits", error.map,ignore.case=TRUE))) {
    #One file provided without .fits extension - must be filelist
    error.map<-try(c(t(read.table(file.path(path.root,error.map), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(error.map)=="try-error") {
      param.warnings<-c(param.warnings,"ErrorMap Parameter table read failed; Using 'NONE'")
      error.map<-"NONE"
    }
    if (is.na(error.map)) {
      param.warnings<-c(param.warnings,"ErrorMap Parameter not in Parameter File; Using 'NONE'")
      error.map<-"NONE"
    }
  }
  #}}}
  #}}}

  #Name of Mask Map {{{
  ID="MaskMap"
  ind<-which(params[ID,]!="")
  mask.map<-params[ID,ind]
  if ((length(ind)==0)||(is.na(mask.map))) {
    param.warnings<-c(param.warnings,"MaskMap Parameter not in Parameter File; Using 'NONE'")
    mask.map<-"NONE"
  }
  #Determine if provided mask.map is an image or filelist {{{
  if ((length(mask.map)==1)&(mask.map!="NONE")&(!grepl(".fits", mask.map,ignore.case=TRUE))) {
    #One file provided without .fits extension - must be filelist
    mask.map<-try(c(t(read.table(file.path(path.root,mask.map), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(mask.map)=="try-error") {
      param.warnings<-c(param.warnings,"MaskMap Parameter table read failed; Using 'NONE'")
      mask.map<-"NONE"
    }
    if (is.na(mask.map)) {
      param.warnings<-c(param.warnings,"MaskMap Parameter not in Parameter File; Using 'NONE'")
      mask.map<-"NONE"
    }
  }
  #}}}
  #}}}

  #Name of Weight Map {{{
  ID="WeightMap"
  ind<-which(params[ID,]!="")
  weight.map<-params[ID,ind]
  if ((length(ind)==0)||(is.na(weight.map))) {
    param.warnings<-c(param.warnings,"WeightMap Parameter not in Parameter File; Using 'NONE'")
    weight.map<-"NONE"
  }
  #Determine if provided weightmap is an image or filelist {{{
  if ((length(weight.map)==1)&(weight.map!="NONE")&(!grepl(".fits", weight.map,ignore.case=TRUE))) {
    #One file provided without .fits extension - must be filelist
    weight.map<-try(c(t(read.table(file.path(path.root,weight.map), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(weight.map)=="try-error") {
      param.warnings<-c(param.warnings,"WeightMap Parameter table read failed; Using 'NONE'")
      weight.map<-"NONE"
    }
    if (is.na(weight.map)) {
      param.warnings<-c(param.warnings,"WeightMap Parameter not in Parameter File; Using 'NONE'")
      weight.map<-"NONE"
    }
  }
  #}}}
  #}}}

  #Zero Point of Weight Map {{{
  ID="WeightMapZP"
  ind<-which(params[ID,]!="")
  wgt.zp<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(wgt.zp))) {
    param.warnings<-c(param.warnings,"WeightMapZP Parameter not in Parameter File; Using 0")
    wgt.zp<-0
  }#}}}

  #Extension number of Data in FITS Header {{{
  ID="DataExtn"
  ind<-which(params[ID,]!="")
  data.extn<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(data.extn))) {
    if ((length(ind)==1)) {
      data.extn<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(data.extn)=='try-error') {
        param.warnings<-c(param.warnings,"DataExtn Parameter table read failed; Using 0")
        data.extn<-0
      }
      if (is.na(data.extn)) {
        param.warnings<-c(param.warnings,"DataExtn Parameter not in Parameter File; Using 0")
        data.extn<-0
      }
    } else {
      param.warnings<-c(param.warnings,"DataExtn Parameter not in Parameter File; Using 0")
      data.extn<-0
    }
  }
  #}}}

  #Extension number of Error Map in FITS Header {{{
  ID="ErrorExtn"
  ind<-which(params[ID,]!="")
  data.error.extn<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(data.error.extn))) {
    if ((length(ind)==1)) {
      data.error.extn<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(data.error.extn)=='try-error') {
        param.warnings<-c(param.warnings,"ErrorExtn Parameter table read failed; Using 0")
        data.error.extn<-0
      }
      if (is.na(data.error.extn)) {
        param.warnings<-c(param.warnings,"ErrorExtn Parameter not in Parameter File; Using 0")
        data.error.extn<-0
      }
    } else {
      param.warnings<-c(param.warnings,"MaskExtn Parameter not in Parameter File; Using 0")
      data.error.extn<-0
    }
  }
  #}}}

  #Extension number of Mask Map in FITS Header {{{
  ID="MaskExtn"
  ind<-which(params[ID,]!="")
  data.mask.extn<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(data.mask.extn))) {
    if ((length(ind)==1)) {
      data.mask.extn<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(data.mask.extn)=='try-error') {
        param.warnings<-c(param.warnings,"MaskExtn Parameter table read failed; Using 0")
        data.mask.extn<-0
      }
      if (is.na(data.mask.extn)) {
        param.warnings<-c(param.warnings,"MaskExtn Parameter not in Parameter File; Using 0")
        data.mask.extn<-0
      }
    } else {
      param.warnings<-c(param.warnings,"MaskExtn Parameter not in Parameter File; Using 0")
      data.mask.extn<-0
    }
  }
  #}}}

  #Extension number of Data in FITS Header {{{
  ID="WeightExtn"
  ind<-which(params[ID,]!="")
  data.weight.extn<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(data.weight.extn))) {
    if ((length(ind)==1)) {
      data.weight.extn<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(data.weight.extn)=='try-error') {
        param.warnings<-c(param.warnings,"WeightExtn Parameter table read failed; Using 0")
        data.weight.extn<-0
      }
      if (is.na(data.weight.extn)) {
        param.warnings<-c(param.warnings,"WeightExtn Parameter not in Parameter File; Using 0")
        data.weight.extn<-0
      }
    } else {
      param.warnings<-c(param.warnings,"WeightExtn Parameter not in Parameter File; Using 0")
      data.weight.extn<-0
    }
  }
  #}}}

  #Do we want to force use of point sources? {{{
  ID="PointSources"
  ind<-which(params[ID,]!="")
  force.point.sources<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(force.point.sources)) {
    if ((length(ind)==1)) {
      force.point.sources<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(force.point.sources)=='try-error') {
        param.warnings<-c(param.warnings,"PointSources Parameter table read failed; Using 0 (FALSE)")
        force.point.sources<-0
      }
      if (is.na(force.point.sources)) {
        param.warnings<-c(param.warnings,"PointSources Parameter not in Parameter File; Using 0 (FALSE)")
        force.point.sources<-0
      }
    } else {
      param.warnings<-c(param.warnings,"PointSources Parameter not in Parameter File; Using 0 (FALSE)")
      force.point.sources<-0
    }

  }
  force.point.sources<-(force.point.sources==1)
  #}}}

  #Error Map scale factor #{{{
  ID="EFactor"
  ind<-which(params[ID,]!="")
  error.factor<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(error.factor))) {
    if ((length(ind)==1)) {
      error.factor<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(error.factor)=='try-error') {
        param.warnings<-c(param.warnings,"EFactor Parameter table read failed; Using 0")
        error.factor<-1
      }
      if (is.na(error.factor)) {
        param.warnings<-c(param.warnings,"EFactor Parameter not in Parameter File; Using 0")
        error.factor<-1
      }
    } else {
      param.warnings<-c(param.warnings,"EFactor Parameter not in Parameter File; Using 0")
      error.factor<-1
    }
  }
  #}}}

  #Flux Correction (Scale) Factor {{{
  ID="FluxCorr"
  ind<-which(params[ID,]!="")
  flux.corr<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(flux.corr))) {
    if ((length(ind)==1)) {
      flux.corr<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(flux.corr)=='try-error') {
        param.warnings<-c(param.warnings,"FluxCorr Parameter table read failed; Using 1")
        flux.corr<-1
      }
      if (is.na(flux.corr)) {
        param.warnings<-c(param.warnings,"FluxCorr Parameter not in Parameter File; Using 1")
        flux.corr<-1
      }
    } else {
      param.warnings<-c(param.warnings,"FluxCorr Parameter not in Parameter File; Using 1")
      flux.corr<-1
    }
  }
  #}}}

  #Initialise Cropping parameters {{{
  crop.radius<-NULL
  ra0<-NULL
  dec0<-NULL
  data.fits.output.filename<-NULL
  mask.fits.output.filename<-NULL
  weight.fits.output.filename<-NULL
  error.fits.output.filename<-NULL
  #}}}

  #Do we want to crop the input image(s)? {{{
  ID="CropImage"
  ind<-which(params[ID,]!="")
  crop.image<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(crop.image))) {
    if ((length(ind)==1)) {
      crop.image<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(crop.image)=='try-error') {
        param.warnings<-c(param.warnings,"CropImage Parameter table read failed; Using 0 (FALSE)")
        crop.image<-0
      }
      if (is.na(crop.image)) {
        param.warnings<-c(param.warnings,"CropImage Parameter not in Parameter File; Using 0 (FALSE)")
        crop.image<-0
      }
    } else {
      param.warnings<-c(param.warnings,"CropImage Parameter not in Parameter File; Using 0 (FALSE)")
      crop.image<-0
    }
  }
  crop.image<-(crop.image==1)
  #}}}

  #Cropped image parameters {{{
  if (crop.image) {
    #What will the cropped image(s) be named {{{
    ID="CropFitsName"
    data.fits.output.filename<-params[ID,1]
    if (is.na(data.fits.output.filename)) {
      param.warnings<-c(param.warnings,"CropFitsName Parameter not in Parameter File; Using 'croppedimage'")
      data.fits.output.filename<-"croppedimage"
    }
    mask.fits.output.filename<-paste(data.fits.output.filename,"_mask.fits",sep="")
    weight.fits.output.filename<-paste(data.fits.output.filename,"_wgt.fits",sep="")
    error.fits.output.filename<-paste(data.fits.output.filename,"_err.fits",sep="")
    data.fits.output.filename<-paste(data.fits.output.filename,".fits",sep="")
    #}}}
    #Cropped image RA {{{
    ID="CropImRA0"
    ind<-which(params[ID,]!="")
    ra0<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(ra0))) {
      if (length(ind)==1) {
        #Try Reading table:
        ra0<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(ra0)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"CropImRA0 Parameter table read failed; Using -999")
          ra0<- -999
        }
        if (is.na(ra0)) {
          #Warn on Error
          param.warnings<-c(param.warnings,"CropImRA0 Parameter not in parameter file; Using -999")
          ra0<- -999
        }
      } else {
        param.warnings<-c(param.warnings,"CropImRA0 Parameter not in Parameter File; Using -999")
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
        dec0<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(dec0)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"CropImDec0 Parameter table read failed; Using -999")
          dec0<- -999
        }
        if (is.na(dec0)) {
          #Warn on Error
          param.warnings<-c(param.warnings,"CropImDec0 Parameter not in parameter file; Using -999")
          dec0<- -999
        }
      } else {
        param.warnings<-c(param.warnings,"CropImDec0 Parameter not in Parameter File; Using -999")
        dec0<- -999
      }
    }
    #}}}
    #Cropped image radius {{{
    ID="CropImRad"
    ind<-which(params[ID,]!="")
    crop.radius<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(crop.radius))) {
      if (length(ind)==1) {
        #Try Reading table:
        crop.radius<-try(as.numeric(c(t(read.table(file.path(path.root,params[ID,1]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(crop.radius)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"CropImRad Parameter table read failed; Using 0.5")
          crop.radius<- 0.5
        }
        if (is.na(crop.radius)) {
          #Warn on Error
          param.warnings<-c(param.warnings,"CropImRad Parameter not in parameter file; Using 0.5")
          crop.radius<- 0.5
        }
      } else {
        param.warnings<-c(param.warnings,"CropImRad Parameter not in Parameter File; Using 0.5")
        crop.radius<-0.5
      }
    }
    #}}}
  }
  #}}}

  #Confusion noise factor (in Units) {{{
  ID="Confusion_units"
  ind<-which(params[ID,]!="")
  conf<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(conf))) {
    if ((length(ind)==1)) {
      conf<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(conf)=='try-error') {
        param.warnings<-c(param.warnings,"Confusion_units Parameter table read failed; Using 0")
        conf<-0
      }
      if (is.na(conf)) {
        param.warnings<-c(param.warnings,"Confusion_units Parameter not in Parameter File; Using 0")
        conf<-0
      }
    } else {
      param.warnings<-c(param.warnings,"Confusion_units Parameter not in Parameter File; Using 0")
      conf<-0
    }
  }
  #}}}

  #Number of Processors available for computations {{{
  ID="nProcessors"
  num.cores<-as.numeric(params[ID,1])
  if (is.na(num.cores)) {
    param.warnings<-c(param.warnings,"nProcessors Parameter not in Parameter File; Using 1")
    num.cores<-1
  }
  #}}}

  #Angular Offset {{{
  #Is there an offset between the Input Catalogue angles and
  #N0E90 Angular Coordinates?
  ID="AngularOffset"
  ind<-which(params[ID,]!="")
  ang.offset<-params[ID,ind]
  if ((length(ind)==0)||is.na(ang.offset)) {
    if ((length(ind)==1)) {
      ang.offset<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(ang.offset)=='try-error') {
        param.warnings<-c(param.warnings,"AngularOffset Parameter table read failed; Using 0 (N0E90)")
        ang.offset<-0
      }
      if (is.na(ang.offset)) {
        param.warnings<-c(param.warnings,"AngularOffset Parameter not in Parameter File; Using 0 (N0E90)")
        ang.offset<-0
      }
    } else {
      param.warnings<-c(param.warnings,"AngularOffset Parameter not in Parameter File; Using 0 (N0E90)")
      ang.offset<-0
    }
  }
  ang.offset<-ang.offset==1
  #}}}

  #Is the map in Jy per Beam? {{{
  ID="MapUnitsPerBeam"
  ind<-which(params[ID,]!="")
  Jybm<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(Jybm))) {
    if ((length(ind)==1)) {
      Jybm<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(Jybm)=='try-error') {
        param.warnings<-c(param.warnings,"MapUnitsPerBeam Parameter table read failed; Using 0 (FALSE)")
        Jybm<-0
      }
      if (is.na(Jybm)) {
        param.warnings<-c(param.warnings,"MapUnitsPerBeam Parameter not in Parameter File; Using 0 (FALSE)")
        Jybm<-0
      }
    } else {
      param.warnings<-c(param.warnings,"MapUnitsPerBeam Parameter not in Parameter File; Using 0 (FALSE)")
      Jybm<-0
    }
  }
  #}}}

  #Resample Apertures {{{
  #Do we want to perform higher precision integrations of
  #Apertures by resampling around the edges?
  ID="ResampleAper"
  ind<-which(params[ID,]!="")
  resample.aperture<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(resample.aperture))) {
    if ((length(ind)==1)) {
      resample.aperture<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(resample.aperture)=="try-error") {
        param.warnings<-c(param.warnings,"ResampleAper Parameter table read failed; Using 1 (TRUE)")
        resample.aperture<-1
      }
      if (is.na(resample.aperture)) {
        param.warnings<-c(param.warnings,"ResampleAper Parameter not in Parameter File; Using 1 (TRUE)")
        resample.aperture<-1
      }
    } else {
      param.warnings<-c(param.warnings,"ResampleAper Parameter not in Parameter File; Using 1 (TRUE)")
      resample.aperture<-1
    }
  }
  resample.aperture<-resample.aperture==1
  #}}}

  #Resample Parameters {{{
  if (any(resample.aperture)) {
    #What resolution do we want to upscale by? {{{
    ID="ResamplingRes"
    ind<-which(params[ID,]!="")
    resample.upres<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(resample.upres))) {
      if ((length(ind)==1)) {
        resample.upres<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(resample.upres)=='try-error') {
          param.warnings<-c(param.warnings,"ResamplingRes Parameter table read failed; Using 3")
          resample.upres<-3
        }
        if (is.na(resample.upres)) {
          param.warnings<-c(param.warnings,"ResamplingRes Parameter not in Parameter File; Using 3")
          resample.upres<-3
        }
      } else {
        param.warnings<-c(param.warnings,"ResamplingRes Parameter not in Parameter File; Using 3")
        resample.upres<-3
      }
    }
    #}}}
    #How many iterations of upscale do we want? {{{
    ID="ResamplingIters"
    ind<-which(params[ID,]!="")
    resample.iterations<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(resample.iterations))) {
      if ((length(ind)==1)) {
        resample.iterations<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(resample.iterations)=='try-error') {
          param.warnings<-c(param.warnings,"ResamplingIters Parameter table read failed; Using 5")
          resample.iterations<-5
        }
        if (is.na(resample.iterations)) {
          param.warnings<-c(param.warnings,"ResamplingIters Parameter not in Parameter File; Using 5")
          resample.iterations<-5
        }
      } else {
        param.warnings<-c(param.warnings,"ResamplingIters Parameter not in Parameter File; Using 5")
        resample.iterations<-5
      }
    }
    #}}}
  } else {
    #If not - set defaults (#iters=0 performs no resampling) {{{
    resample.upres<-2
    resample.iterations<-0
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
      confidence<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(confidence)=='try-error') {
        param.warnings<-c(param.warnings,"PSFConfidence Parameter table read failed; Using 1")
        confidence<-1
      }
      if (is.na(confidence)) {
        param.warnings<-c(param.warnings,"PSFConfidence Parameter not in Parameter File; Using 1")
        confidence<-1
      }
    } else {
      param.warnings<-c(param.warnings,"PSFConfidence Parameter not in Parameter File; Using 1")
      confidence<-1
    }
  }
  #}}}

  #Size of the aperture stamp as a multiple of the aperture major axis {{{
  ID="ApStampWidth"
  ind<-which(params[ID,]!="")
  def.buff<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(def.buff)) {
    if ((length(ind)==1)) {
      def.buff<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(def.buff)=='try-error') {
        param.warnings<-c(param.warnings,"ApStampWidth Parameter table read failed; Using 1.05")
        def.buff<-1.05
      }
      if (is.na(def.buff)) {
        param.warnings<-c(param.warnings,"ApStampWidth Parameter not in Parameter File; Using 1.05")
        def.buff<-1.05
      }
    } else {
      param.warnings<-c(param.warnings,"ApStampWidth Parameter not in Parameter File; Using 1.05")
      def.buff<-1.05
    }
  }
  if (def.buff<1) {
    param.warnings<-c(param.warnings,"ApStampWidth Value is less than or equal to Unity. Value must be strictly > 1. Setting to 1.05")
    def.buff<-1.05
  }
  #}}}

  #Do we want to output the source mask only? {{{
  ID="SourceMaskOnly"
  ind<-which(params[ID,]!="")
  sourcemask.only<-as.numeric(params[ID,ind])
  if (length(ind)==0||is.na(sourcemask.only)) {
    if (length(ind)==1) {
      sourcemask.only<-try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
      if (class(sourcemask.only)=="try-error") {
        param.warnings<-c(param.warnings,"SourceMaskOnly Parameter table read failed; Using 0 (FALSE)")
        sourcemask.only<-0
      }
      if (is.na(sourcemask.only)) {
        param.warnings<-c(param.warnings,"SourceMaskOnly Parameter not in Parameter File; Using 0 (FALSE)")
        sourcemask.only<-0
      }
    } else {
      param.warnings<-c(param.warnings,"SourceMaskOnly Parameter not in Parameter File; Using 0 (FALSE)")
      sourcemask.only<-0
    }
  }
  sourcemask.only<-(sourcemask.only==1)
  #}}}

  #Do we want to output the source mask at all? {{{
  if (any(!sourcemask.only)) {
    ID="WriteSourceMask"
    ind<-which(params[ID,]!="")
    sourcemask.out<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(sourcemask.out)) {
      if (length(ind)==1) {
        sourcemask.out<-try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
        if (class(sourcemask.out)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"WriteSourceMask Parameter table read failed; Using 0 (FALSE)")
          sourcemask.out<-0
        }
        if (is.na(sourcemask.out)) {
          param.warnings<-c(param.warnings,"WriteSourceMask Parameter not in Parameter File; Using 0 (FALSE)")
          sourcemask.out<-0
        }
      } else {
        param.warnings<-c(param.warnings,"WriteSourceMask Parameter not in Parameter File; Using 0 (FALSE)")
        sourcemask.out<-0
      }
    }
  } else {
    sourcemask.out<-TRUE
    sourcemask   <-TRUE
  }
  sourcemask.out<-(sourcemask.out==1)
  sourcemask   <-(sourcemask.out==1)
  #}}}

  #Do we want to output the All Apertures Mask {{{
  all.apertures.map.filename<-NULL
  ID="WriteAAMask"
  ind<-which(params[ID,]!="")
  make.all.apertures.map<-params[ID,ind]
  if ((length(ind)==0)||is.na(make.all.apertures.map)) {
    param.warnings<-c(param.warnings,"WriteAAMask Parameter not in Parameter File; Using 0 (FALSE)")
    make.all.apertures.map<-FALSE
  } else { make.all.apertures.map<-(make.all.apertures.map==1) }
  #}}}

  #All Apertures mask file name {{{
  if ( make.all.apertures.map ) {
    #Name of the output All Apertures file
    ID="AllApersFile"
    ind<-which(params[ID,]!="")
    all.apertures.map.filename<-params[ID,ind]
    if ((length(ind)==0)||is.na(all.apertures.map.filename)) {
      param.warnings<-c(param.warnings,"AllApersFile Parameter not in Parameter File; Using 'AllApertures_Mask.fits'")
      all.apertures.map.filename<-"AllApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Convolved Apertures mask {{{
  fa.filename<-NULL
  ID="WriteFAMask"
  ind<-which(params[ID,]!="")
  make.convolved.apertures.map<-params[ID,ind]
  if ((length(ind)==0)||is.na(make.convolved.apertures.map)) {
    param.warnings<-c(param.warnings,"WriteFAMask Parameter not in Parameter File; Using 0 (FALSE)")
    make.convolved.apertures.map<-FALSE
  } else { make.convolved.apertures.map<-(make.convolved.apertures.map==1) }
  #}}}

  #Convolved Apertures Filename {{{
  if ( make.convolved.apertures.map ) {
    #Name of the output Convolved Apertures Mask
    ID="ConvApersFile"
    ind<-which(params[ID,]!="")
    fa.filename<-params[ID,ind]
    if ((length(ind)==0)||is.na(fa.filename)) {
      param.warnings<-c(param.warnings,"ConvApersFile Parameter not in Parameter File; Using 'AllConvolvedApertures_Mask.fits'")
      fa.filename<-"AllConvolvedApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Deblended Convolved Apertures mask {{{
  dfa.filename<-NULL
  ID="WriteDFAMask"
  ind<-which(params[ID,]!="")
  make.debelended.apertures.map<-params[ID,ind]
  if ((length(ind)==0)||is.na(make.debelended.apertures.map)) {
    param.warnings<-c(param.warnings,"WriteDFAMask Parameter not in Parameter File; Using 0 (FALSE)")
    make.debelended.apertures.map<-FALSE
  } else { make.debelended.apertures.map<-(make.debelended.apertures.map==1) }
  #}}}

  #Deblended Convolved Apertures Mask filename {{{
  if ( make.debelended.apertures.map ) {
    #Name of the output Convolved Apertures Mask
    ID="DeblConvApersFile"
    ind<-which(params[ID,]!="")
    dfa.filename<-params[ID,ind]
    if ((length(ind)==0)||is.na(dfa.filename)) {
      param.warnings<-c(param.warnings,"DeblConvApersFile Parameter not in Parameter File; Using 'AllDeblConvolvedApertures_Mask.fits'")
      dfa.filename<-"AllDeblConvolvedApertures_Mask.fits"
    }
  }
  #}}}

  #Do we want to output the Residual image? {{{
  residual.map<-NULL
  ID="WriteResidMap"
  ind<-which(params[ID,]!="")
  make.resid.map<-params[ID,ind]
  if ((length(ind)==0)||is.na(make.resid.map)) {
    param.warnings<-c(param.warnings,"WriteResidMap Parameter not in Parameter File; Using 0 (FALSE)")
    make.resid.map<-FALSE
  } else { make.resid.map<-(make.resid.map==1) }
  #}}}

  #Residual Image filename {{{
  if ( make.resid.map ) {
    #Name of the output Residual image
    ID="ResidImageFile"
    ind<-which(params[ID,]!="")
    residual.map<-params[ID,ind]
    if ((length(ind)==0)||is.na(residual.map)) {
      param.warnings<-c(param.warnings,"ResidImageFile Parameter not in Parameter File; Using 'ResidualImage.fits'")
      residual.map<-"ResidualImage.fits"
    }
  }
  #}}}

  #Do we want to output the Flux table? {{{
  tableout.name<-NULL
  ID="WriteTable"
  ind<-which(params[ID,]!="")
  write.tab<-params[ID,ind]
  if ((length(ind)==0)||is.na(write.tab)) {
    param.warnings<-c(param.warnings,"WriteTable Parameter not in Parameter File; Using 1 (TRUE)")
    write.tab<-TRUE
  } else { write.tab<-(write.tab==1) }
  #}}}

  #Table Filename {{{
  if ( write.tab ) {
    #Name of output Flux table
    ID="TableName"
    ind<-which(params[ID,]!="")
    tableout.name<-params[ID,ind]
    if ((length(ind)==0)||(is.na(tableout.name))) {
      #Warn on Error
      param.warnings<-c(param.warnings,"TableName Parameter not in Parameter File; Using 'LAMBDAR_Results'")
      tableout.name<-"LAMBDAR_Results"
    } else {
      if (length(ind)==1) {
        tableout.name<-try(suppressWarnings(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
        if (class(tableout.name)=="try-error") {
          #Warn on Error
          tableout.name<-params[ID,1]
        }
        if (is.na(tableout.name)) {
          tableout.name<-params[ID,1]
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
      param.warnings<-c(param.warnings,"ShowTime Parameter not in Parameter File; Using 0 (FALSE)")
      showtime<-FALSE
    } else { showtime<-(showtime==1) }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Interactive"
  ind<-which(params[ID,]!="")
  interact<-params[ID,ind]
  if ((length(ind)==0)||is.na(interact)) {
    param.warnings<-c(param.warnings,"Interactive Parameter not in Parameter File; Using 0 (FALSE)")
    interact<-FALSE
  } else { interact<-(interact==1) }
  #}}}

  #What limit do we want for the use of masks what cross the mask edges {{{
  ID="UseMaskLim"
  ind<-which(params[ID,]!="")
  use.mask.lim<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(use.mask.lim)) {
    if ((length(ind)==1)) {
      use.mask.lim<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(use.mask.lim)=='try-error') {
        param.warnings<-c(param.warnings,"UseMaskLim Parameter table read failed; Using 0.2")
        use.mask.lim<-0.2
      }
      if (is.na(use.mask.lim)) {
        param.warnings<-c(param.warnings,"UseMaskLim Parameter not in Parameter File; Using 0.2")
        use.mask.lim<-0.2
      }
    } else {
      param.warnings<-c(param.warnings,"UseMaskLim Parameter not in Parameter File; Using 0.2")
      use.mask.lim<-0.2
    }
  }
  #}}}

  #Do we want Diagnostic Output in Log File {{{
  ID="Diagnostic"
  ind<-which(params[ID,]!="")
  diagnostic<-params[ID,ind]
  if ((length(ind)==0)||is.na(diagnostic)) {
    param.warnings<-c(param.warnings,"Diagnostic Parameter not in Parameter File; Using 0 (FALSE)")
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
      param.warnings<-c(param.warnings,"Verbose Parameter not in Parameter File; Using 0 (FALSE)")
      verbose<-FALSE
    } else { verbose<-(verbose==1) }
  }
  verbose.out<-verbose
  if (quiet) {
    verbose<-FALSE
    diagnostic<-FALSE
  }
  #}}}

  #Do we want a sample of the apertures to be output? {{{
  ID="PlotSample"
  ind<-which(params[ID,]!="")
  plot.sample<-params[ID,ind]
  if ((length(ind)==0)||is.na(plot.sample)) {
    param.warnings<-c(param.warnings,"PlotSample Parameter not in Parameter File; Using 0 (FALSE)")
    plot.sample<-FALSE
  } else { plot.sample<-(plot.sample==1) }
  ID="PlotAll"
  ind<-which(params[ID,]!="")
  plot.all<-params[ID,ind]
  if ((length(ind)==0)||is.na(plot.all)) {
    param.warnings<-c(param.warnings,"PlotAll Parameter not in Parameter File; Using 0 (FALSE)")
    plot.all<-FALSE
  } else { plot.all<-(plot.all==1) }
  if (plot.all & !plot.sample) {
    param.warnings<-c(param.warnings,"PlotAll Parameter TRUE but PlotSample Parameter FALSE; Forcing PlotSample Parameter to TRUE")
    plot.sample<-TRUE
  }
  #}}}

  #Make magnitudes in Output? {{{
  ID="Magnitudes"
  ind<-which(params[ID,]!="")
  magnitudes<-params[ID,ind]
  if ((length(ind)==0)||is.na(magnitudes)) {
    param.warnings<-c(param.warnings,"Magnitudes Parameter not present in the Parameter File; Using 1 (TRUE)")
    magnitudes<-TRUE
  } else { magnitudes<-(magnitudes==1) }
  #}}}

  #Magnitude Details {{{
  if (magnitudes) {
    #AB Vega Magnitude {{{
    ID="ABVegaFlux"
    ind<-which(params[ID,]!="")
    ab.vega.flux<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(ab.vega.flux))) {
      param.warnings<-c(param.warnings,"ABVegaFlux Parameter not present in the Parameter File; Using 1.0")
      ab.vega.flux<-1.0
    }
    #}}}

    #magnitudes Zero Point {{{
    ID="MagZeroPoint"
    ind<-which(params[ID,]!="")
    mag.zp<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(mag.zp))) {
      if (length(ind)==1) {
        mag.zp<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(mag.zp)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"MagZeroPoint Parameter table read failed; Using 0.0")
          mag.zp<-0.0
        }
        if (is.na(mag.zp)) {
          #Warn on Error
          param.warnings<-c(param.warnings,"MagZeroPoint Parameter not present/bad in the Parameter File; Using 0.0")
          mag.zp<-0.0
        }
      } else {
        #Warn on Error
        param.warnings<-c(param.warnings,"MagZeroPoint Parameter not present/bad in the Parameter File; Using 0.0")
        mag.zp<-0.0
      }
    }
    #}}}

    #magnitudes Zero Point {{{
    ID="MagZPLabel"
    ind<-which(params[ID,]!="")
    mag.zp.label<-params[ID,ind]
    if ((length(ind)==0)||(is.na(mag.zp.label))) {
      param.warnings<-c(param.warnings,"MagZPLabel Parameter not present in the Parameter File; Using 'MagZP'")
      mag.zp.label<-"MagZP"
    }
    #}}}
  } else {
    mag.zp<-NA
    mag.zp.label<-"MagZP"
    ab.vega.flux<-1.0
  }
  #}}}

  #Saturation Label {{{
  ID="SaturationLabel"
  ind<-which(params[ID,]!="")
  satur.label<-params[ID,ind]
  if ((length(ind)==0)||(is.na(satur.label))) {
    param.warnings<-c(param.warnings,"SaturationLabel Parameter not present in the Parameter File; Using 'SATUR'")
    satur.label<-"SATUR"
  } else if (length(ind)==1&&(file.exists(file.path(path.root,satur.label)))) {
    #One file provided
    satur.label<-try(c(t(read.table(file.path(path.root,satur.label), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(satur.label)=="try-error") {
      param.warnings<-c(param.warnings,"SaturationLabel Parameter table read failed; Using 'SATUR'")
      satur.label<-"SATUR"
    }
    if (is.na(satur.label)) {
      param.warnings<-c(param.warnings,"SaturationLabel Parameter not in Parameter File; Using 'SATUR'")
      satur.label<-"SATUR"
    }
  }
  #}}}

  #Saturation Label {{{
  ID="Saturation"
  ind<-which(params[ID,]!="")
  saturation<-as.numeric(params[ID,ind])
  if ((length(ind)==0)) {
    param.warnings<-c(param.warnings,"Saturation Parameter not present in the Parameter File; Using Inf")
    saturation<-Inf
  } else if (length(ind)==1&&(file.exists(file.path(path.root,saturation)))) {
    #One file provided
    saturation<-try(as.numeric(c(t(read.table(file.path(path.root,saturation), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#")))),silent=TRUE)
    if (class(saturation)=="try-error") {
      param.warnings<-c(param.warnings,"Saturation Parameter table read failed; Using Inf")
      saturation<-Inf
    }
    if (is.na(saturation)) {
      param.warnings<-c(param.warnings,"Saturation Parameter not in Parameter File; Using Inf")
      saturation<-Inf
    }
  } else if (any(is.na(saturation))) {
    param.warnings<-c(param.warnings,"Saturation Parameter is NA in Parameter File; Using Inf")
    saturation[which(is.na(saturation))]<-Inf
  }
  #}}}

  #Gain Label {{{
  ID="GainLabel"
  ind<-which(params[ID,]!="")
  gain.label<-params[ID,ind]
  if ((length(ind)==0)||(is.na(gain.label))) {
    param.warnings<-c(param.warnings,"GainLabel Parameter not present in the Parameter File; Using 'GAIN'")
    gain.label<-"GAIN"
  } else if (length(ind)==1&&(file.exists(file.path(path.root,gain.label)))) {
    #One file provided
    gain.label<-try(c(t(read.table(file.path(path.root,gain.label), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE)
    if (class(gain.label)=="try-error") {
      param.warnings<-c(param.warnings,"GainLabel Parameter table read failed; Using 'GAIN'")
      gain.label<-"GAIN"
    }
    if (is.na(gain.label)) {
      param.warnings<-c(param.warnings,"GainLabel Parameter not in Parameter File; Using 'GAIN'")
      gain.label<-"GAIN"
    }
  }
  #}}}

  #Blanks Correction {{{
  #Do we want to perform a blanks Correction to the errors of object fluxes?
  ID="BlankCor"
  ind<-which(params[ID,]!="")
  blank.cor<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(blank.cor))) {
    if ((length(ind)==1)) {
      blank.cor<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(blank.cor)=="try-error") {
        param.warnings<-c(param.warnings,"BlankCor Parameter table read failed; Using 0 (FALSE)")
        blank.cor<-0
      }
      if (is.na(blank.cor)) {
        param.warnings<-c(param.warnings,"BlankCor Parameter not in Parameter File; Using 0 (FALSE)")
        blank.cor<-0
      }
    } else {
      param.warnings<-c(param.warnings,"BlankCor Parameter not in Parameter File; Using 0 (FALSE)")
      blank.cor<-0
    }
  }
  blank.cor<-blank.cor==1
  #Number of Blanks per Object {{{
  ID="nBlanks"
  ind<-which(params[ID,]!="")
  num.blanks<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(num.blanks))) {
    if ((length(ind)==1)) {
      num.blanks<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(num.blanks)=="try-error") {
        param.warnings<-c(param.warnings,"nBlanks Parameter table read failed; Using 10")
        num.blanks<-10
      }
      if (is.na(num.blanks)) {
        param.warnings<-c(param.warnings,"nBlanks Parameter not in Parameter File; Using 10")
        num.blanks<-10
      }
    } else {
      param.warnings<-c(param.warnings,"nBlanks Parameter not in Parameter File; Using 10")
      num.blanks<-10
    }
  }
  #}}}
  #}}}

  #Randoms Correction {{{
  #Do we want to perform a randoms Correction to the errors of object fluxes?
  ID="RanCor"
  ind<-which(params[ID,]!="")
  ran.cor<-params[ID,ind]
  if (length(ind)!=0){ if(ran.cor!="execute") { ran.cor<-as.numeric(ran.cor) } }
  if ((length(ind)==0)||(is.na(ran.cor))) {
    if ((length(ind)==1)) {
      ran.cor<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(ran.cor)=="try-error") {
        param.warnings<-c(param.warnings,"RanCor Parameter table read failed; Using 0 (FALSE)")
        ran.cor<-0
      }
      if (is.na(ran.cor)) {
        param.warnings<-c(param.warnings,"RanCor Parameter not in Parameter File; Using 0 (FALSE)")
        ran.cor<-0
      }
    } else {
      param.warnings<-c(param.warnings,"RanCor Parameter not in Parameter File; Using 0 (FALSE)")
      ran.cor<-0
    }
  }
  if(ran.cor!="execute") { ran.cor<-ran.cor==1 }
  #Number of Randoms per Object {{{
  ID="nRandoms"
  ind<-which(params[ID,]!="")
  num.randoms<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(num.randoms))) {
    if ((length(ind)==1)) {
      num.randoms<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(num.randoms)=="try-error") {
        param.warnings<-c(param.warnings,"nRandoms Parameter table read failed; Using 10")
        num.randoms<-10
      }
      if (is.na(num.randoms)) {
        param.warnings<-c(param.warnings,"nRandoms Parameter not in Parameter File; Using 10")
        num.randoms<-10
      }
    } else {
      param.warnings<-c(param.warnings,"nRandoms Parameter not in Parameter File; Using 10")
      num.randoms<-10
    }
  }
  #}}}
  #}}}

  #Perform a Sky estimation & subtraction? {{{
  ID="DoSkyEst"
  ind<-which(params[ID,]!="")
  do.sky.est<-as.numeric(params[ID,ind])
  if ((length(ind)==0) || is.na(do.sky.est)) {
    if (length(ind)==1) {
      do.sky.est<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(do.sky.est)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"DoSkyEst Parameter table read failed; Using 0 (FALSE)")
        do.sky.est<-0
      }
      if (is.na(do.sky.est)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"DoSkyEst Parameter not in Parameter File; Using 0 (FALSE)")
        do.sky.est<-0
      }
    } else {
      #Warn on Error
      param.warnings<-c(param.warnings,"DoSkyEst Parameter not in Parameter File; Using 0 (FALSE)")
      do.sky.est<-0
    }
  }
  do.sky.est<-(do.sky.est==1)
  #}}}

  #Calculate the Sky RMS? {{{
  ID="GetSkyRMS"
  ind<-which(params[ID,]!="")
  get.sky.rms<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(get.sky.rms)) {
    if (length(ind)==1) {
      get.sky.rms<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(get.sky.rms)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"GetSkyRMS Parameter not in Parameter File; Using 0 (FALSE)")
        get.sky.rms<-0
      }
      if (is.na(get.sky.rms)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"GetSkyRMS Parameter not in Parameter File; Using 0 (FALSE)")
        get.sky.rms<-0
      }
    } else {
        #Warn on Error
        param.warnings<-c(param.warnings,"GetSkyRMS Parameter not in Parameter File; Using 0 (FALSE)")
        get.sky.rms<-0
    }
  }
  get.sky.rms<-get.sky.rms==1
  #}}}

  #Sky Estimate Paramters {{{
  if (any(do.sky.est|get.sky.rms|blank.cor)) {
    #Sourcemask needed for SkyEstimate. If not TRUE, set to TRUE {{{
    if (any(!sourcemask & (do.sky.est|get.sky.rms|blank.cor))) {
      param.warnings<-c(param.warnings,"Source Mask creation being forced for all runs with doSkyEst/getSkyRMS/blank.cor Parameters TRUE")
      sourcemask<-blank.cor|do.sky.est|get.sky.rms|sourcemask
    }
    #}}}

    #Number of iterations used in sky estimation {{{
    ID="SkyEstIters"
    ind<-which(params[ID,]!="")
    sky.clip.iters<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||is.na(sky.clip.iters)) {
      param.warnings<-c(param.warnings,"SkyEstIters Parameter not in Parameter File; Using 5")
      sky.clip.iters<-5
    }
    #}}}

    #Sigma level used in sky cut {{{
    ID="SkyEstProbCut"
    ind<-which(params[ID,]!="")
    sky.clip.prob<-as.numeric(params[ID,ind])
    if ((length(ind)==0)|| (is.na(sky.clip.prob))) {
      if ((length(ind)==1)) {
        sky.clip.prob<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(sky.clip.prob)=='try-error') {
          param.warnings<-c(param.warnings,"SkyEstProbCut Parameter table read failed; Using 3")
          sky.clip.prob<-3
        }
        if (is.na(sky.clip.prob)) {
          param.warnings<-c(param.warnings,"SkyEstProbCut Parameter not in Parameter File; Using 3")
          sky.clip.prob<-3
        }
      } else {
        param.warnings<-c(param.warnings,"SkyEstProbCut Parameter not in Parameter File; Using 3")
        sky.clip.prob<-3
      }
    }
    #}}}

    #Default Sky Value if estimation fails {{{
    ID="SkyDefault"
    ind<-which(params[ID,]!="")
    sky.default<-params[ID,ind]
    if ((length(ind)==0)||is.na(sky.default)) {
      param.warnings<-c(param.warnings,"SkyDefault Parameter not in Parameter File; Using 'median'")
      sky.default<-"median"
    }
    #}}}

    #Level of Corellated noise in image {{{
    ID="SkyCorrelNoise"
    ind<-which(params[ID,]!="")
    correl.noise<-as.numeric(params[ID,ind])
    if ((length(ind)==0)||(is.na(correl.noise))) {
      if ((length(ind)==1)) {
        correl.noise<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(correl.noise)=='try-error') {
          param.warnings<-c(param.warnings,"SkyCorrelNoise Parameter table read failed; Using 1")
          correl.noise<-1
        }
        if (is.na(correl.noise)) {
          param.warnings<-c(param.warnings,"SkyCorrelNoise Parameter not in Parameter File; Using 1")
          correl.noise<-1
        }
      } else {
        param.warnings<-c(param.warnings,"SkyCorrelNoise Parameter not in Parameter File; Using 1")
        correl.noise<-1
      }
    }
    #}}}
  } else {
    sky.clip.iters<-0
    sky.default<-0
    sky.clip.prob<-0
    correl.noise<-1
  }
  #}}}

  #Sourcemask parameters; filename & confidence limit {{{
  sourcemask.filename<-NULL
  sourcemask.conf.lim<-NULL
  if ( sourcemask ) {
    if (sourcemask.out) {
      #Name of SourceMask that is output
      ID="SourceMaskFile"
      ind<-which(params[ID,]!="")
      sourcemask.filename<-params[ID,ind]
      if (length(ind)==0||is.na(sourcemask.filename)||((length(sourcemask.filename)==1)&(!grepl(".fits", sourcemask.filename,ignore.case=TRUE)))) {
        if (length(ind)==1) {
          sourcemask.filename<-(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
          if (class(sourcemask.filename)=="try-error") {
            #Warn on Error
            param.warnings<-c(param.warnings,"SourceMaskFile Parameter table read failed; Using 'SourceMask.fits'")
            sourcemask.filename<-"SourceMask.fits"
          }
          if (is.na(sourcemask.filename)) {
            #Warn on Error
            param.warnings<-c(param.warnings,"SourceMaskFile Parameter not in Parameter File; Using 'SourceMask.fits'")
            sourcemask.filename<-"SourceMask.fits"
          }
        } else {
          #Warn on Error
        param.warnings<-c(param.warnings,"SourceMaskFile Parameter not in Parameter File; Using 'SourceMask.fits'")
          sourcemask.filename<-"SourceMask.fits"
        }
      }
    }
    #Name of SourceMask that is output
    ID="TransmissionMap"
    ind<-which(params[ID,]!="")
    transmission.map<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(transmission.map)) {
      if (length(ind)==1) {
        transmission.map<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(transmission.map)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"TransmissionMap Parameter table read failed; Using 0 (FALSE)")
          transmission.map<-0
        }
        if (is.na(transmission.map)) {
          #Warn on Error
          param.warnings<-c(param.warnings,"TransmissionMap Parameter not in Parameter File; Using 0 (FALSE)")
          transmission.map<-0
        }
      } else {
        #Warn on Error
        param.warnings<-c(param.warnings,"TransmissionMap Parameter not in Parameter File; Using 0 (FALSE)")
        transmission.map<-0
      }
    }
    transmission.map<-(transmission.map==1)
    #SourceMask Confidence Limit
    ID="SourceMaskConfLim"
    ind<-which(params[ID,]!="")
    sourcemask.conf.lim<-as.numeric(params[ID,ind])
    if (length(ind)==0||is.na(sourcemask.conf.lim)) {
      if (length(ind)==1) {
        sourcemask.conf.lim<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
        if (class(sourcemask.conf.lim)=="try-error") {
          #Warn on Error
          param.warnings<-c(param.warnings,"SourceMaskConfLim Parameter table read failed; Using 0.95")
          sourcemask.conf.lim<-0.95
        }
        if (is.na(sourcemask.conf.lim)) {
          #Warn on Error
          param.warnings<-c(param.warnings,"SourceMaskConfLim Parameter not in Parameter File; Using 0.95")
          sourcemask.conf.lim<-0.95
        }
      } else {
        #Warn on Error
        param.warnings<-c(param.warnings,"SourceMaskConfLim Parameter not in Parameter File; Using 0.95")
        sourcemask.conf.lim<-0.95
      }
    }
  } else {
    transmission.map=FALSE
    sourcemask.conf.lim=NA
  }
  #}}}

  #Set Minimum Aperture Radius {{{
  ID="MinApRad"
  ind<-which(params[ID,]!="")
  min.ap.rad<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(min.ap.rad))) {
    if (length(ind)==1) {
      min.ap.rad<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(min.ap.rad)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"MinApRad Parameter not in Parameter File; Using 0")
        min.ap.rad<-0
      }
      if (is.na(min.ap.rad)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"MinApRad Parameter not in Parameter File; Using 0")
        min.ap.rad<-0
      }
    } else {
      #Warn on Error
      param.warnings<-c(param.warnings,"MinApRad Parameter not in Parameter File; Using 0")
      min.ap.rad<-0
    }
  }
  #}}}

  #Do we want a Memory-Safe run? {{{
  ID="MemorySafe"
  ind<-which(params[ID,]!="")
  mem.safe<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(mem.safe)) {
    if ((length(ind)==1)) {
      mem.safe<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(mem.safe)=='try-error') {
        param.warnings<-c(param.warnings,"MemorySafe Parameter table read failed; Using 0 (FALSE)")
        mem.safe<-0
      }
      if (is.na(mem.safe)) {
        mem.safe<-0
        param.warnings<-c(param.warnings,"MemorySafe Parameter not in Parameter File; Using 0 (FALSE)")
      }
    } else {
      mem.safe<-0
      param.warnings<-c(param.warnings,"MemorySafe Parameter not in Parameter File; Using 0 (FALSE)")
    }
  }
  if (mem.safe > 1) { force.safe=TRUE } else { force.safe=FALSE }
  mem.safe<-mem.safe >= 1
  #}}}

  #What limit value do we want for the Aperture generation? #{{{
  ID="ApertureConfLimit"
  ind<-which(params[ID,]!="")
  ap.limit<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||(is.na(ap.limit))) {
    if (length(ind)==1) {
      ap.limit<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(ap.limit)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"ApertureConfLimit Parameter table read failed; Using 0.9")
        ap.limit<-0.9
      }
      if (is.na(ap.limit)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"ApertureConfLimit Parameter not in Parameter File; Using 0.9")
        ap.limit<-0.9
      }
    } else {
      #Warn on Error
      param.warnings<-c(param.warnings,"ApertureConfLimit Parameter not in Parameter File; Using 0.9")
      ap.limit<-0.9
    }
  }
  #}}}

  #Do we want to iterate the fluxes? {{{
  ID="IterateFluxes"
  ind<-which(params[ID,]!="")
  iterate.fluxes<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(iterate.fluxes)) {
    if (length(ind)==1) {
      iterate.fluxes<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(iterate.fluxes)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"IterateFluxes Parameter table read failed; Using 0 (FALSE)")
        iterate.fluxes<-0
      }
      if (is.na(iterate.fluxes)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"IterateFluxes Parameter not in Parameter File; Using 0 (FALSE)")
        iterate.fluxes<-0
      }
    } else {
        #Warn on Error
        param.warnings<-c(param.warnings,"IterateFluxes Parameter not in Parameter File; Using 0 (FALSE)")
        iterate.fluxes<-0
    }
  }
  iterate.fluxes<-(iterate.fluxes==1)
  #}}}

  #If Iterating, how many iterations do we want? {{{
  ID="nIterations"
  ind<-which(params[ID,]!="")
  num.iterations<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(num.iterations)) {
    if ((length(ind)==1)) {
      num.iterations<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(num.iterations)=='try-error') {
        param.warnings<-c(param.warnings,"nIterations Parameter table read failed; Using 2")
        num.iterations<-2
      }
      if (is.na(num.iterations)) {
        param.warnings<-c(param.warnings,"nIterations Parameter not in Parameter File; Using 2")
        num.iterations<-2
      }
    } else {
      param.warnings<-c(param.warnings,"nIterations Parameter not in Parameter File; Using 2")
      num.iterations<-2
    }
  }
  #}}}

  #What format are the input flux.weights? {{{
  ID="FluxWgtType"
  ind<-which(params[ID,]!="")
  weight.type<-tolower(params[ID,ind])
  if ((length(ind)==0)||is.na(weight.type)) {
    param.warnings<-c(param.warnings,"FluxWgtType Parameter not in Parameter File; Using 'scale'")
    weight.type<-"scale"
  } else if (weight.type!="flux" & weight.type!="mag" & weight.type!="scale") {
    stop("Fluxweight Type Parameter has unknown value; valid values are 'flux', 'mag', and 'scale'")
  }
  #}}}

  #Do we want to flux.weight Using Pixel Fluxes? {{{
  ID="PixelFluxWgt"
  ind<-which(params[ID,]!="")
  use.pixel.fluxweight<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(use.pixel.fluxweight)) {
    if (length(ind)==1) {
      use.pixel.fluxweight<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(use.pixel.fluxweight)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"PixelFluxWgt Parameter table read failed; Using 0 (FALSE)")
        use.pixel.fluxweight<-0
      }
      if (is.na(use.pixel.fluxweight)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"PixelFluxWgt Parameter not in Parameter File; Using 0 (FALSE)")
        use.pixel.fluxweight<-0
      }
    } else {
      #Warn on Error
      param.warnings<-c(param.warnings,"PixelFluxWgt Parameter not in Parameter File; Using 0 (FALSE)")
      use.pixel.fluxweight<-0
    }
  }
  use.pixel.fluxweight<-(use.pixel.fluxweight==1)
  if (use.pixel.fluxweight) { weight.type<-"scale" }
  #}}}

  #Do we just want the Deblend Fraction for each aperture? {{{
  ID="GetDeblFrac"
  ind<-which(params[ID,]!="")
  get.debl.frac<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(get.debl.frac)) {
    if (length(ind)==1) {
      get.debl.frac<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(get.debl.frac)=="try-error") {
        #Warn on Error
        param.warnings<-c(param.warnings,"GetDeblFrac Parameter not in Parameter File; Using 0 (FALSE)")
        get.debl.frac<-0
      }
      if (is.na(get.debl.frac)) {
        #Warn on Error
        param.warnings<-c(param.warnings,"GetDeblFrac Parameter not in Parameter File; Using 0 (FALSE)")
        get.debl.frac<-0
      }
    } else {
        #Warn on Error
        param.warnings<-c(param.warnings,"GetDeblFrac Parameter not in Parameter File; Using 0 (FALSE)")
        get.debl.frac<-0
    }
  }
  get.debl.frac<-get.debl.frac==1
  #Don't do unwanted things if using this option {{{
  if (any(get.debl.frac & (get.sky.rms|do.sky.est))) {
    param.warnings<-c(param.warnings,"Stopping Sky Estimate for all loops where GetDeblFrac is TRUE, as it is unnecessary")
    if (length(do.sky.est)==1) {
      do.sky.est<-!get.debl.frac
    } else if (length(get.debl.frac)==1) {
      do.sky.est<-FALSE
    } else {
      do.sky.est[which(get.debl.frac)]<-FALSE
    }
    if (length(get.sky.rms)==1) {
      get.sky.rms<-!get.debl.frac
    } else if (length(get.debl.frac)==1) {
      get.sky.rms<-FALSE
    } else {
      get.sky.rms[which(get.debl.frac)]<-FALSE
    }
  }
  #}}}
  #}}}

  #Simulation Noise Gaussianisation Kernel {{{
  ID="SimGauss_AS"
  ind<-which(params[ID,]!="")
  sim.gauss.arcsec<-as.numeric(params[ID,ind])
  if ((length(ind)==0)||is.na(sim.gauss.arcsec)) {
    if ((length(ind)==1)) {
      sim.gauss.arcsec<-as.numeric(try(c(t(read.table(file.path(path.root,params[ID,ind[1]]), strip.white=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE, comment.char = "#"))),silent=TRUE))
      if (class(sim.gauss.arcsec)=='try-error') {
        param.warnings<-c(param.warnings,"SimGauss_AS Parameter table read failed; Using 0")
        sim.gauss.arcsec<-0
      }
      if (is.na(num.iterations)) {
        param.warnings<-c(param.warnings,"SimGauss_AS Parameter not in Parameter File; Using 0")
        sim.gauss.arcsec<-0
      }
    } else {
      param.warnings<-c(param.warnings,"SimGauss_AS Parameter not in Parameter File; Using 0")
      sim.gauss.arcsec<-0
    }
  }
  #}}}

  #Name of Logfile to be output {{{
  ID="LogFile"
  logfile<-params[ID,1]
  if (is.na(logfile)) {
    param.warnings<-c(param.warnings,"LogFile Parameter not in Parameter File; Using 'LAMBDAR_Log.txt'")
    logfile<-"LAMBDAR_Log.txt"
  }
  #}}}
  #}}}

  # Print any warnings {{{
  if (!is.null(param.warnings) & !quiet) {
    param.warnings<-paste(param.warnings,collapse="\n     > ")
    cat("{\n    Warnings in Parameter File read:\n     > ")
    cat(param.warnings)
    cat("\n   } ")
  }

  #}}}

  # Assign variables to LAMBDAR workspace {{{
  assign("all.apertures.map.filename"       , all.apertures.map.filename       , envir = env) # A
  assign("ab.vega.flux"       , ab.vega.flux       , envir = env) #
  assign("ang.offset"        , ang.offset        , envir = env) #
  assign("ap.limit"          , ap.limit          , envir = env) #
  assign("beam.area.input.as"  , beam.area.input.as  , envir = env) # B
  assign("blank.cor"         , blank.cor         , envir = env) # B
  assign("conf"             , conf             , envir = env) # C
  assign("check.contam"      , check.contam      , envir = env) #
  assign("confidence"       , confidence       , envir = env) #
  assign("crop.image"        , crop.image        , envir = env) #
  assign("correl.noise"     , correl.noise     , envir = env) #
  assign("catalogue"        , catalogue        , envir = env) #
  assign("crop.radius"           , crop.radius           , envir = env) #
  assign("contam.lab"        , contam.lab        , envir = env) #
  assign("cata.lab"          , cata.lab          , envir = env) #
  assign("data.map"          , data.map          , envir = env) # D
  assign("do.sky.est"         , do.sky.est         , envir = env) #
  assign("def.buff"          , def.buff          , envir = env) #
  assign("diagnostic"       , diagnostic       , envir = env) #
  assign("dec0"             , dec0             , envir = env) #
  assign("dec.lab"           , dec.lab           , envir = env) #
  assign("dfa.filename"      , dfa.filename      , envir = env) #
  assign("error.map"         , error.map         , envir = env) # E
  assign("data.extn"             , data.extn             , envir = env) #
  assign("data.error.extn"          , data.error.extn          , envir = env) #
  assign("data.mask.extn"         , data.mask.extn         , envir = env) #
  assign("data.weight.extn"          , data.weight.extn          , envir = env) #
  assign("error.factor"          , error.factor          , envir = env) #
  assign("flux.corr"         , flux.corr         , envir = env) # F
  assign("fa.filename"       , fa.filename       , envir = env) #
  assign("force.point.sources", force.point.sources, envir = env) #
  assign("force.safe", force.safe, envir = env) #
  assign("filt.contam"       , filt.contam       , envir = env) #
  assign("flux.weight.lab"    , flux.weight.lab    , envir = env) #
  assign("gain.label"        , gain.label        , envir = env) #
  assign("gauss.fwhm.arcsec"    , gauss.fwhm.arcsec    , envir = env) # G
  assign("get.debl.frac"      , get.debl.frac      , envir = env) #
  assign("get.sky.rms"        , get.sky.rms        , envir = env) #
  assign("group.weights"       , group.weights       , envir = env) #
  assign("resample.iterations"        , resample.iterations        , envir = env) # I
  assign("iterate.fluxes"    , iterate.fluxes    , envir = env) # I
  assign("num.iterations"      , num.iterations      , envir = env) # I
  assign("interact"         , interact         , envir = env) #
  assign("mask.fits.output.filename"   , mask.fits.output.filename   , envir = env) #
  assign("weight.fits.output.filename" , weight.fits.output.filename , envir = env) #
  assign("error.fits.output.filename"   , error.fits.output.filename   , envir = env) #
  assign("data.fits.output.filename"    , data.fits.output.filename    , envir = env) #
  assign("Jybm"             , Jybm             , envir = env) # J
  assign("logfile"          , logfile          , envir = env) # KL
  assign("make.resid.map"     , make.resid.map     , envir = env) # M
  assign("make.debelended.apertures.map"      , make.debelended.apertures.map      , envir = env) #
  assign("magnitudes"       , magnitudes       , envir = env) #
  assign("mag.zp"            , mag.zp            , envir = env) #
  assign("mag.zp.label"       , mag.zp.label       , envir = env) #
  assign("make.convolved.apertures.map"       , make.convolved.apertures.map       , envir = env) #
  assign("make.all.apertures.map"       , make.all.apertures.map       , envir = env) #
  assign("mask.map"          , mask.map          , envir = env) #
  assign("min.ap.rad"         , min.ap.rad         , envir = env) #
  assign("mem.safe"          , mem.safe          , envir = env) #
  assign("no.psf"            , no.psf            , envir = env) # N
  assign("num.nearest.neighbours"             , num.nearest.neighbours             , envir = env) #
  assign("no.contam.map"      , no.contam.map      , envir = env) #
  assign("num.cores"           , num.cores           , envir = env) #
  assign("num.blanks"          , num.blanks          , envir = env) #
  assign("num.randoms"         , num.randoms         , envir = env) #
  assign("optimal.aper"         , optimal.aper         , envir = env) # O
  assign("path.root"         , path.root         , envir = env) # P
  assign("path.work"         , path.work         , envir = env) #
  assign("path.out"          , path.out          , envir = env) #
  assign("plot.sample"       , plot.sample       , envir = env) #
  assign("plot.all"          , plot.all          , envir = env) #
  assign("psf.map"           , psf.map           , envir = env) #
  assign("psf.weighted"      , psf.weighted      , envir = env) #
  assign("psf.filt"          , psf.filt          , envir = env) #
  assign("psf.label"         , psf.label         , envir = env) #
  assign("psf.label.type"    , psf.label.type    , envir = env) #
  assign("resample.aperture" , resample.aperture , envir = env) # QR
  assign("ra0"              , ra0              , envir = env) #
  assign("ra.lab"            , ra.lab            , envir = env) #
  assign("ran.cor"           , ran.cor           , envir = env) #
  assign("residual.map"         , residual.map         , envir = env) #
  assign("saturation"       , saturation       , envir = env) # S
  assign("satur.label"       , satur.label       , envir = env) #
  assign("sim.gauss.arcsec"     , sim.gauss.arcsec     , envir = env) #
  assign("sourcemask"       , sourcemask       , envir = env) #
  assign("sourcemask.out"    , sourcemask.out    , envir = env) #
  assign("sourcemask.only"   , sourcemask.only   , envir = env) #
  assign("sourcemask.conf.lim"  , sourcemask.conf.lim  , envir = env) #
  assign("showtime"         , showtime         , envir = env) #
  assign("sky.clip.iters"      , sky.clip.iters      , envir = env) #
  assign("sky.default"       , sky.default       , envir = env) #
  assign("sky.clip.prob"       , sky.clip.prob       , envir = env) #
  assign("semimaj.lab"       , semimaj.lab       , envir = env) #
  assign("semimin.lab"       , semimin.lab       , envir = env) #
  assign("sourcemask.filename"       , sourcemask.filename       , envir = env) #
  assign("tableout.name"     , tableout.name     , envir = env) # T
  assign("theta.lab"         , theta.lab         , envir = env) #
  assign("transmission.map"  , transmission.map  , envir = env) #
  assign("resample.upres"            , resample.upres            , envir = env) # U
  assign("use.mask.lim"       , use.mask.lim       , envir = env)
  assign("use.pixel.fluxweight", use.pixel.fluxweight , envir = env)
  assign("verbose"          , verbose          , envir = env) # V
  assign("verbose.out"       , verbose.out       , envir = env) #
  assign("write.tab"         , write.tab         , envir = env) # W
  assign("weight.type"       , weight.type       , envir = env) #
  assign("weight.map"           , weight.map           , envir = env) #
  assign("wgt.zp"            , wgt.zp            , envir = env) #
                                                              # XYZ
  #}}}

  #Remove unneeded variables {{{
  rm(params)
  #}}}

  #Finished Setup of Parameter Space {{{
  if (!quiet) { cat(" - Done\n") }
  return=param.warnings
  #}}}

}

opencatalogue <-
function(outenv=parent.env(environment()), env=NULL){

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

  #Open and Read Catalogue {{{
  if (!quiet) { cat(paste("   Reading Input Catalogue",catalogue,"   ")) }
  if (showtime) { timer<-proc.time() }

  #Determine catalogue type {{{
  csv=FALSE
  fits=FALSE
  rdat=FALSE
  if (grepl(".csv", catalogue)) { csv=TRUE }
  else if (grepl(".fits", catalogue)) { fits=TRUE }
  else if (grepl(".Rdata", catalogue)) { rdat=TRUE }
  else { stop("Catalogue does not have a recognised extension (.csv/.fits/.Rdata)") }
  #}}}
  #Open Catalogue {{{
  if (csv) {
    fitstable<-try(as.data.frame(read.csv(paste(pathroot,pathwork,catalogue,sep=""),stringsAsFactors=F)))
    #Test Read of Catalogue for errors
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (fits) {
    fitstable<-try(as.data.frame(read.fitstab(paste(pathroot,pathwork,catalogue,sep="")),stringsAsFactors=F))
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (rdat) {
    names<-try(as.list(load(paste(pathroot,pathwork,catalogue,sep=""))))
    if (class(names)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
    fitstable<-get(names[[1]])
  }
  #}}}
  #Get catalogue size {{{
  nrows<-length(fitstable[,1])
  #}}}
  #}}}

  #Check for Correct Column Syntax & Read Data {{{
  #Catalogue ID {{{
  id_g<-try(fitstable[1:nrows,catalab])
  if ((class(id_g)=="try-error")||(is.null(id_g[1]))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",catalab,"column"))
  }#}}}
  #Make sure that there are no duplicate IDs {{{
  if (any(is.na(id_g))|any(id_g=="",na.rm=T)) {
    message("There are",length(which(is.na(id_g)|id_g=="")),"Missing IDs. These will be renamed NewID_%d")
    ind<-which(is.na(id_g)|id_g=="")
    id_g[ind]<-paste("NewID_",cbind(1:length(ind)),sep="")
  }
  if (any(duplicated(id_g))) {
    message("There are",length(which(duplicated(id_g))),"duplicated IDs. These will be appended with DuplicID_%d")
    ind<-which(duplicated(id_g))
    id_g[ind]<-paste(id_g[ind],"_DuplicID_",c(1:length(ind)),sep="")
  }
  #}}}
  #Object RA {{{
  ra_g<-try(as.numeric(fitstable[1:nrows,ralab]))
  if ((class(ra_g)=="try-error")||(is.null(ra_g[1]))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",ralab,"column"))
  }#}}}
  #Object Dec {{{
  dec_g<-try(as.numeric(fitstable[1:nrows,declab]))
  if ((class(dec_g)=="try-error")||(is.null(dec_g[1]))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",declab,"column"))
  }#}}}
  #Setup Aperture Variables {{{
  if (forcepointsources) {
    #If forcing point sources, set standardised aperture values {{{
    theta_g<-array(0, dim=c(nrows))
    a_g<-array(0, dim=c(nrows))
    b_g<-array(0, dim=c(nrows))
    #}}}
  } else {
    #Otherwise, Check Syntax & Read in Aperture Variables {{{
    #Aperture Angle {{{
    theta_g<-try(as.numeric(fitstable[1:nrows,thetalab]))  # theta
    if ((class(theta_g)=="try-error")||(is.null(theta_g))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",thetalab,"column"))
    }#}}}
    #Aperture Semi-Major Axis {{{
    a_g<-try(as.numeric(fitstable[1:nrows,semimajlab])) # semimajor in arcsec
    if ((class(a_g)=="try-error")||(is.null(a_g))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",semimajlab,"column"))
    }
    a_g[which(!is.finite(a_g))]<-0
    #}}}
    #Aperture Semi-Minor Axis {{{
    b_g<-try(as.numeric(fitstable[1:nrows,semiminlab])) # semiminor in arcsec
    if ((class(b_g)=="try-error")||(is.null(b_g))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",semiminlab,"column"))
    }
    b_g[which(!is.finite(b_g))]<-0
    #}}}
    #}}}
  }
  #}}}

  #If wanted, read contaminants column {{{
  if (filtcontam) {
    contams<-try(as.numeric(fitstable[1:nrows,contamlab]))
    if ((class(contams)=="try-error")||(is.null(contams))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",contamlab,"column"))
    }
    message(paste("There are",length(which(contams==1)),"Contaminants to be subtracted"))
  }#}}}

  #If Weight Column exists, read values {{{
  fluxweight<-try(as.numeric(fitstable[1:nrows,fluxweightlab] ))
  if ((class(fluxweight)=="try-error")||(is.null(fluxweight))) {
    #Otherwise, set all weights to unity
    fluxweight<-1
  }#}}}
  #}}}

  #Parse Parameter Space {{{
  assign("id_g"      ,id_g      ,envir=outenv)
  assign("ra_g"      ,ra_g      ,envir=outenv)
  assign("dec_g"     ,dec_g     ,envir=outenv)
  assign("a_g"       ,a_g       ,envir=outenv)
  assign("b_g"       ,b_g       ,envir=outenv)
  assign("theta_g"   ,theta_g   ,envir=outenv)
  if (filtcontam) { assign("contams"   ,contams   ,envir=outenv) }
  assign("fluxweight",fluxweight,envir=outenv)
  #assign("fitstable" ,fitstable ,envir=outenv)
  assign("nrows"     ,nrows     ,envir=outenv)
  #}}}

  #Finished Reading Catalogue, return {{{
  if (showtime) { cat(paste(" - Done (", round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  } else if (!quiet) { cat(" - Done\n") }
  if (!is.null(env)) { detach(env) }
  return=NULL
  #}}}

}

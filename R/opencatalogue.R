opencatalogue <-
function(env=NULL, outenv=NULL){

  # Load Parameter Space
  if (is.null(env)) {
    stop("No Parameter Space Environment Specified in Call")
  }
  if (is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))

  #Open and Read Catalogue
  if (!quiet) { cat(paste("   Reading Input Catalogue",catalogue,"   ")) }
  if (showtime) { timer<-proc.time() }

  csv=FALSE
  fits=FALSE
  rdat=FALSE
  if (grepl(".csv", catalogue)) { csv=TRUE }
  else if (grepl(".fits", catalogue)) { fits=TRUE }
  else if (grepl(".Rdata", catalogue)) { rdat=TRUE }
  else { stop("Catalogue does not have a recognised extension (.csv/.fits/.Rdata)") }
  if (csv) {
    error<-try(as.data.frame(read.csv(paste(pathroot,catalogue,sep=""))))
    #Test Read of Catalogue for errors
    if (class(error)=="try error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
    fitstable<-as.data.frame(read.csv(paste(pathroot,catalogue,sep="")))
  } else if (fits) {
    error<-try(as.data.frame(read.fitstab(paste(pathroot,catalogue,sep=""))))
    if (class(error)=="try error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
    fitstable<-as.data.frame(read.fitstab(paste(pathroot,catalogue,sep="")))
  } else if (rdat) {
    error<-try(load(paste(pathroot,catalogue,sep="")))
    if (class(error)=="try error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
    names<-as.list(load(paste(pathroot,catalogue,sep="")))
    fitstable<-get(names[[1]])
  }
  #Read Fits Table
  nrows<-length(fitstable[,1])
  #Check for Correct Column Syntax & Read Data
  #GAMA Catalogue ID
  error<-try(id_g<-fitstable[1:nrows,"CATAID"])
  if ((class(error)=="try-error")||(is.null(id_g[1]))) {
    sink(type="message")
    stop("Catalogue does not contain CATAID column")
  }
  #Object RA
  error<-try(ra_g<-fitstable[1:nrows,"RA"])
  if ((class(error)=="try-error")||(is.null(ra_g[1]))) {
    sink(type="message")
    stop("Catalogue does not contain RA column")
  }
  #Object Dec
  error<-try(dec_g<-fitstable[1:nrows,"DEC"])
  if ((class(error)=="try-error")||(is.null(dec_g[1]))) {
    sink(type="message")
    stop("Catalogue does not contain DEC column")
  }

  #Setup Aperture Variables
  if (forcepointsources) {
    #If using point - set standardised aperture values
    theta_g<-array(0, dim=c(nrows))
    a_g<-array(0, dim=c(nrows))
    b_g<-array(0, dim=c(nrows))
  } else {
    #Otherwise, Check Syntax & Read in Aperture Variables
    error<-try(theta_g<-fitstable[1:nrows,"THETA_J2000"])  # theta
    if ((class(error)=="try-error")||(is.null(theta_g))) {
      error<-try(theta_g<-fitstable[1:nrows,"THETA_J2000_N0E90"])  # theta
      if ((class(error)=="try-error")||(is.null(theta_g))) { sink(type="message") ; stop("Catalogue does not contain THETA_J2000 column") }
    }
    error<-try(a_g<-fitstable[1:nrows,"A_KRON_SEMIMAJ_ARC"]) # semimajor in arcsec
    if ((class(error)=="try-error")||(is.null(a_g))) { sink(type="message") ; stop("Catalogue does not contain A_KRON_SEMIMAJ_ARC column") }
    a_g[which(!is.finite(a_g))]<-0
    error<-try(b_g<-fitstable[1:nrows,"B_KRON_SEMIMIN_ARC"]) # semiminor in arcsec
    if ((class(error)=="try-error")||(is.null(b_g))) { sink(type="message") ; stop("Catalogue does not contain B_KRON_SEMIMIN_ARC column") }
    b_g[which(!is.finite(b_g))]<-0
  }

  #Are we removing contaminants?
  if (filtcontam) {
    #If so, check that contaminant flag column exists
    error<-try(contams<-fitstable[1:nrows,"CONTAM"] )
    if ((class(error)=="try-error")||(is.null(contams))) { sink(type="message") ; stop("Catalogue does not contain CONTAM column") }
    message(paste("There are",length(which(contams==1)),"Contaminants to be subtracted"))
  }

  #If Weight Column exists, read values
  error<-try(fluxweight<-fitstable[1:nrows,"WEIGHT"] )
  if ((class(error)=="try-error")||(is.null(fluxweight))) {
    #Otherwise, set all weights to unity
    fluxweight<-1
  }

  #Parse Parameter Space
  assign("id_g"      ,id_g      ,envir=outenv)
  assign("ra_g"      ,ra_g      ,envir=outenv)
  assign("dec_g"     ,dec_g     ,envir=outenv)
  assign("a_g"       ,a_g       ,envir=outenv)
  assign("b_g"       ,b_g       ,envir=outenv)
  assign("theta_g"   ,theta_g   ,envir=outenv)
  if (filtcontam) { assign("contams"   ,contams   ,envir=outenv) }
  assign("fluxweight",fluxweight,envir=outenv)
  assign("fitstable" ,fitstable ,envir=outenv)
  assign("nrows"     ,nrows     ,envir=outenv)


  #Finished Reading Catalogue
  if (showtime) { cat(paste(" - Done (", round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  } else if (!quiet) { cat(" - Done\n") }
  detach(env)

}

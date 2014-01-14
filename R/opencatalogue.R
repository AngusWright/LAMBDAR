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
    fitstable<-try(as.data.frame(read.csv(paste(pathroot,catalogue,sep=""))))
    #Test Read of Catalogue for errors
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (fits) {
    fitstable<-try(as.data.frame(read.fitstab(paste(pathroot,catalogue,sep=""))))
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (rdat) {
    names<-try(as.list(load(paste(pathroot,catalogue,sep=""))))
    if (class(names)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
    fitstable<-get(names[[1]])
  }
  #Read Fits Table
  nrows<-length(fitstable[,1])
  #Check for Correct Column Syntax & Read Data
  #GAMA Catalogue ID
  id_g<-try(fitstable[1:nrows,catalab])
  if ((class(id_g)=="try-error")||(is.null(id_g[1]))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",catalab,"column"))
  }
  #Object RA
  ra_g<-try(fitstable[1:nrows,ralab])
  if ((class(ra_g)=="try-error")||(is.null(ra_g[1]))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",ralab,"column"))
  }
  #Object Dec
  dec_g<-try(fitstable[1:nrows,declab])
  if ((class(dec_g)=="try-error")||(is.null(dec_g[1]))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",declab,"column"))
  }

  #Setup Aperture Variables
  if (forcepointsources) {
    #If using point - set standardised aperture values
    theta_g<-array(0, dim=c(nrows))
    a_g<-array(0, dim=c(nrows))
    b_g<-array(0, dim=c(nrows))
  } else {
    #Otherwise, Check Syntax & Read in Aperture Variables
    #Aperture Angle
    theta_g<-try(fitstable[1:nrows,thetalab])  # theta
    if ((class(theta_g)=="try-error")||(is.null(theta_g))) {
      sink(type="message") 
      stop(paste("Catalogue does not contain",thetalab,"column"))
    }
    #Aperture Semi-Major Axis
    a_g<-try(fitstable[1:nrows,semimajlab]) # semimajor in arcsec
    if ((class(a_g)=="try-error")||(is.null(a_g))) { 
      sink(type="message") 
      stop(paste("Catalogue does not contain",semimajlab,"column"))
    }
    a_g[which(!is.finite(a_g))]<-0
    #Aperture Semi-Minor Axis
    b_g<-try(fitstable[1:nrows,semiminlab]) # semiminor in arcsec
    if ((class(b_g)=="try-error")||(is.null(b_g))) { 
      sink(type="message") 
      stop(paste("Catalogue does not contain",semiminlab,"column"))
    }
    b_g[which(!is.finite(b_g))]<-0
  }

  #Are we removing contaminants?
  if (filtcontam) {
    #If so, check that contaminant flag column exists
    contams<-try(fitstable[1:nrows,contamlab] )
    if ((class(contams)=="try-error")||(is.null(contams))) { 
      sink(type="message") 
      stop(paste("Catalogue does not contain",contamlab,"column")) 
    }
    message(paste("There are",length(which(contams==1)),"Contaminants to be subtracted"))
  }

  #If Weight Column exists, read values
  fluxweight<-try(fitstable[1:nrows,fluxweightlab] )
  if ((class(fluxweight)=="try-error")||(is.null(fluxweight))) {
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

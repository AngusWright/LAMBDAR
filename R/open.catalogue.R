open.catalogue <-
function(outenv=parent.env(environment()), save.table=FALSE, env=NULL){

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
  ascii=FALSE
  fits=FALSE
  rdat=FALSE
  if (grepl(".csv", catalogue,ignore.case=TRUE)) { csv=TRUE }
  else if (grepl(".dat", catalogue,ignore.case=TRUE)) { ascii=TRUE }
  else if (grepl(".fits", catalogue,ignore.case=TRUE)) { fits=TRUE }
  else if (grepl(".Rdata", catalogue,ignore.case=TRUE)) { rdat=TRUE }
  else { stop("Catalogue does not have a recognised extension (.csv/.dat/.fits/.Rdata)") }
  #}}}
  #Open Catalogue {{{
  if (csv) {
    fitstable<-try(fread(paste(path.root,path.work,catalogue,sep=""),data.table=FALSE,stringsAsFactors=FALSE),silent=TRUE)
    #Test Read of Catalogue for errors
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (ascii) {
    fitstable<-try(as.data.frame(read.table(paste(path.root,path.work,catalogue,sep="")),stringsAsFactors=FALSE),silent=TRUE)
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (fits) {
    fitstable<-try(read.fits.cat(paste(path.root,path.work,catalogue,sep=""),data.table=FALSE,stringsAsFactors=FALSE),silent=TRUE)
    if (class(fitstable)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
  } else if (rdat) {
    names<-try(as.list(load(paste(path.root,path.work,catalogue,sep=""))),silent=TRUE)
    if (class(names)=="try-error") {
      #Stop on Error
      geterrmessage()
      sink(type="message")
      stop("Catalogue File read failed")
    }
    fitstable<-get(names[[1]])
    if (!is.data.frame(fitstable)) { 
      fitstable<-as.data.frame(fitstable)
    }
  }
  #}}}
  #Get catalogue size {{{
  num.rows<-length(fitstable[,1])
  #}}}
  #}}}

  #Check for Correct Column Syntax & Read Data {{{
  #Catalogue ID {{{
  cat.id<-try(fitstable[,cata.lab],silent=TRUE)
  if ((class(cat.id)=="try-error")||(length(cat.id)==0)||all(is.na(cat.id))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",cata.lab,"column"))
  }#}}}
  #Make sure that there are no duplicate IDs {{{
  if (any(is.na(cat.id))|any(cat.id=="",na.rm=TRUE)) {
    message("There are ",length(which(is.na(cat.id)|cat.id==""))," Missing IDs. These will be renamed NewID_%d")
    ind<-which(is.na(cat.id)|cat.id=="")
    cat.id[ind]<-paste("NewID_",cbind(1:length(ind)),sep="")
  }
  if (any(duplicated(cat.id))) {
    message("There are ",length(which(duplicated(cat.id)))," duplicated IDs. These will be appended with DuplicID_%d")
    ind<-which(duplicated(cat.id))
    cat.id[ind]<-paste(cat.id[ind],"_DuplicID_",c(1:length(ind)),sep="")
  }
  #}}}
  #Object RA {{{
  cat.ra<-try(as.numeric(fitstable[,ra.lab]),silent=TRUE)
  if ((class(cat.ra)=="try-error")||(length(cat.ra)==0)||all(is.na(cat.ra))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",ra.lab,"column"))
  }#}}}
  #Object Dec {{{
  cat.dec<-try(as.numeric(fitstable[,dec.lab]),silent=TRUE)
  if ((class(cat.dec)=="try-error")||(length(cat.dec)==0)||all(is.na(cat.dec))) {
    sink(type="message")
    stop(paste("Catalogue does not contain",dec.lab,"column"))
  }#}}}
  #Setup Aperture Variables {{{
  if (force.point.sources) {
    #If forcing point sources, set standardised aperture values {{{
    cat.theta<-array(0, dim=c(num.rows))
    cat.a<-array(0, dim=c(num.rows))
    cat.b<-array(0, dim=c(num.rows))
    #}}}
  } else {
    #Otherwise, Check Syntax & Read in Aperture Variables {{{
    #Aperture Angle {{{
    cat.theta<-try(as.numeric(fitstable[,theta.lab]),silent=TRUE)  # theta
    if ((class(cat.theta)=="try-error")||(length(cat.theta)==0)||all(is.na(cat.theta))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",theta.lab,"column"))
    }#}}}
    #Aperture Semi-Major Axis {{{
    cat.a<-try(as.numeric(fitstable[,semimaj.lab]),silent=TRUE) # semimajor in arcsec
    if ((class(cat.a)=="try-error")||(length(cat.a)==0)||all(is.na(cat.a))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",semimaj.lab,"column"))
    }
    cat.a[which(!is.finite(cat.a))]<-0
    #}}}
    #Aperture Semi-Minor Axis {{{
    cat.b<-try(as.numeric(fitstable[,semimin.lab]),silent=TRUE) # semiminor in arcsec
    if ((class(cat.b)=="try-error")||(length(cat.b)==0)||all(is.na(cat.b))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",semimin.lab,"column"))
    }
    cat.b[which(!is.finite(cat.b))]<-0
    #}}}
    #}}}
  }
  #}}}
  #If save.table, recreate the catalogue with the checked-parameters {{{
  if (save.table) { 
    fitstable[,cata.lab]<-cat.id
    fitstable[,ra.lab]<-cat.ra
    fitstable[,dec.lab]<-cat.dec
    fitstable[,semimaj.lab]<-cat.a
    fitstable[,semimin.lab]<-cat.b
    fitstable[,theta.lab]<-cat.theta
  }
  #}}}
  #If wanted, read contaminants column {{{
  if (filt.contam) {
    contams<-try(as.numeric(fitstable[,contam.lab]),silent=TRUE)
    if ((class(contams)=="try-error")||(length(contams)==0)||all(is.na(contams))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",contam.lab,"column"))
    }
    message(paste("There are ",length(which(contams==1))," Contaminants to be subtracted"))
    #If save.table, recreate the catalogue with the checked-parameters {{{
    if (save.table) { 
      fitstable[,contam.lab]<-contams
    }
    #}}}
  } else { 
    #If save.table, recreate the catalogue with the checked-parameters {{{
    if (save.table) { 
      fitstable[,contam.lab]<-rep(0,num.rows)
    }
    #}}}
  }#}}}
  #If Weight Column exists, read values {{{
  flux.weight<-try(as.numeric(fitstable[,flux.weight.lab] ),silent=TRUE)
  if ((class(flux.weight)=="try-error")||(length(flux.weight)==0)||all(is.na(flux.weight))) {
    #Otherwise, set all weights to unity
    flux.weight<-1
  }#}}}

  #If wanted, read grouping column {{{
  if (!exists("group.weights")) { group.weights<-FALSE }
  if (group.weights) {
    if (!exists("group.lab")) { group.lab<-"GROUP" }
    groups<-try(as.numeric(fitstable[1:num.rows,group.lab]),silent=TRUE)
    if ((class(groups)=="try-error")||(length(groups)==0)||all(is.na(groups))) {
      sink(type="message")
      stop(paste("Catalogue does not contain",group.lab,"column"))
    }
    message(paste("There are ",length(factor(groups))," Groups being used in weighting"))
  }#}}}
  #If save.table, recreate the catalogue with the checked-parameters {{{
  if (save.table) { 
    fitstable[,flux.weight.lab]<-flux.weight
    #Remove unneeded columns
    fitstable<-fitstable[,c(cata.lab,ra.lab,dec.lab,semimaj.lab,semimin.lab,theta.lab,flux.weight.lab,contam.lab)]
  }
  #}}}
  #}}}

  #Parse Parameter Space {{{
  assign("cat.id"      ,cat.id      ,envir=outenv)
  assign("cat.ra"      ,cat.ra      ,envir=outenv)
  assign("cat.dec"     ,cat.dec     ,envir=outenv)
  assign("cat.a"       ,cat.a       ,envir=outenv)
  assign("cat.b"       ,cat.b       ,envir=outenv)
  assign("cat.theta"   ,cat.theta   ,envir=outenv)
  if (filt.contam) { assign("contams"   ,contams   ,envir=outenv) }
  if (group.weights) { assign("groups"   ,groups   ,envir=outenv) }
  assign("flux.weight",flux.weight,envir=outenv)
  assign("num.rows"     ,num.rows     ,envir=outenv)
  if (save.table) { 
    assign("saved.table"     ,fitstable     ,envir=outenv)
  }
  #}}}

  #Finished Reading Catalogue, return {{{
  if (showtime) { cat(paste(" - Done (", round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  } else if (!quiet) { cat(" - Done\n") }
  if (!is.null(env)) { detach(env) }
  return=NULL
  #}}}

}

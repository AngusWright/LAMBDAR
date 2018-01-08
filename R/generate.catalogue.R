generate.catalogue <-
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

  #Open and Read Segmentation image {{{
  if (!quiet) { cat(paste("   Reading Segmentation Image",catalogue,"   ")) }
  if (showtime) { timer<-proc.time() }

  #Determine catalogue type {{{
  if (grepl(".fits", catalogue,ignore.case=TRUE)) { fits=TRUE }
  else { stop("Segmentation Map does not have a recognised extension (.fits)") }
  #}}}
  #Read the Image astrometry {{{
  seg.astr<-try(read.astrometry(paste(path.root,path.work,catalogue,sep="")),silent=TRUE)
  if (class(seg.astr)=="try-error") {
    #Stop on Error
    geterrmessage()
    sink(type="message")
    stop("Segmentation Image read failed")
  }
  #}}}
  #Read the Image {{{
  segmentation<-try(read.fits.im(paste(path.root,path.work,catalogue,sep=""))$dat[[1]],silent=TRUE)
  if (class(segmentation)=="try-error") {
    #Stop on Error
    geterrmessage()
    sink(type="message")
    stop("Segmentation Image read failed")
  }
  #}}}
  #Catalogue ID {{{
  cat.id<-sort(unique(as.numeric(segmentation)))[-1]
  #}}}
  #Get catalogue size {{{
  num.rows<-length(cat.id)
  #}}}
  #}}}

  #Check for Correct Column Syntax & Read Data {{{
  #Object RA & DEC {{{
  seg.pix<-abs(seg.astr$CD[1,1])*3600
  cat.xy<-foreach(i=cat.id,.combine='rbind',.inorder=TRUE,
                  .export=c("segmentation",'colMedians','colMaxs','colMins'))%dopar%{
    inds=which(segmentation==i,arr.ind=TRUE)
    return=c(colMedians(inds),(colMaxs(inds)-colMins(inds))/2)
  }
  seg.x<-cat.xy[,1]
  seg.y<-cat.xy[,2]
  cat.pos<-xy.to.ad(seg.x,seg.y,astr.struc=seg.astr)
  cat.ra<-cat.pos[,1]
  cat.dec<-cat.pos[,2]
  cat.a<-ifelse(cat.xy[,3]>cat.xy[,4],cat.xy[,3],cat.xy[,4])*seg.pix
  cat.b<-ifelse(cat.xy[,3]>cat.xy[,4],cat.xy[,4],cat.xy[,3])*seg.pix
  cat.theta<-ifelse(cat.xy[,3]>cat.xy[,4],0,90)
  #}}}
  #If save.table, recreate the catalogue with the checked-parameters {{{
  if (save.table) {
    fitstable<-data.frame(cat.id,cat.ra,cat.dec,cat.a,cat.b,cat.theta,seg.x,seg.y)
    segx.lab<-"SegIm_X"
    segy.lab<-"SegIm_Y"
    colnames(fitstable)<-c(cata.lab,ra.lab,dec.lab,semimaj.lab,semimin.lab,theta.lab,segx.lab,segy.lab)
  }
  #}}}
  #Initialise contaminants column {{{
  contams<-rep(0,num.rows)
  message(paste("Apertures from Segmentation map, there are no contaminants!"))
  #If save.table, recreate the catalogue with the checked-parameters {{{
  if (save.table) {
    fitstable[,contam.lab]<-contams
  }
  #}}}
  #}}}
  #Initialise flux-weights column {{{
  flux.weight<-1
  if (save.table) {
    fitstable[,flux.weight.lab]<-flux.weight
  }
  #}}}

  #If wanted, read grouping column {{{
  if (!exists("group.weights")) { group.weights<-FALSE }
  if (!exists("group.lab")) { group.lab<-"GROUP" }
  message(paste("Using segmentation map, there can be no Groups used in weighting"))
  groups<-rep(1,num.rows)
  #}}}
  #If save.table, recreate the catalogue with the checked-parameters {{{
  if (save.table) {
    fitstable[,group.lab]<-groups
  }
  #}}}
  #}}}

  #Parse Parameter Space {{{
  assign("cat.id"      ,cat.id      ,envir=outenv)
  assign("cat.ra"      ,cat.ra      ,envir=outenv)
  assign("cat.dec"     ,cat.dec     ,envir=outenv)
  assign("seg.x"       ,seg.x       ,envir=outenv)
  assign("seg.y"       ,seg.y       ,envir=outenv)
  assign("cat.a"       ,cat.a       ,envir=outenv)
  assign("cat.b"       ,cat.b       ,envir=outenv)
  assign("cat.theta"   ,cat.theta   ,envir=outenv)
  if (filt.contam) { assign("contams"   ,contams   ,envir=outenv) }
  if (group.weights) { assign("groups"   ,groups   ,envir=outenv) }
  assign("flux.weight",flux.weight,envir=outenv)
  assign("num.rows"     ,num.rows     ,envir=outenv)
  assign("segmentation" ,segmentation    ,envir=outenv)
  assign("seg.astr" ,seg.astr    ,envir=outenv)
  assign("seg.pix" ,seg.pix    ,envir=outenv)
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

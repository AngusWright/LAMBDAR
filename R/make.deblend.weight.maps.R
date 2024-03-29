make.deblend.weight.maps <-
function(outenv=parent.env(environment()), env=NULL,subs=NULL) {
#Procedure produces an array of stamps with apertures that
#have been weighted according to their amount of blending.
#This provides a launchpad for the final determination of
#per-object fluxes.

  if (!quiet) { cat("Make_Deblended_Weightmap    ") }
  message('--------------------Make_Deblended_Weightmap---------------------------')

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

  #Setup Sizes {{{
  if (is.null(subs)) {
    ngal<-length(cat.id)
    subs<-1:ngal
  } else {
    ngal<-length(subs)
  }
  #}}}

  #Produce Deblended Stamps {{{
  #Details {{{
  #For each stamp, determine the limits of the stamp
  #in full (weighted) mask array-space, and divide stamp
  #by the full value. }}}
  #Perform Calculation and assignment {{{
  if (cutup) { registerDoParallel(cores=1) } 
  deblend.weight.map<-foreach(wsfam=wsfa[subs], i=1:ngal,xlo=ap.lims.data.map[subs,1],xup=ap.lims.data.map[subs,2],ylo=ap.lims.data.map[subs,3],yup=ap.lims.data.map[subs,4], .export="image.env$wfa", .inorder=TRUE)%dopar%{
      #Check for errors and create deblended stamps {{{
      if ((sum(wsfam) > 0)&&(sum(image.env$wfa[xlo:xup,ylo:yup])>0)) {
        #Aperture stamps are error free {{{
        #Check for errors in array sizes {{{
        if (length(wsfam[1,]) != length(image.env$wfa[xlo:xup,ylo:yup][1,])) {
          print(paste(dim(wsfam),"!=",dim(image.env$wfa[xlo:xup,ylo:yup])))
          print(paste("[",xlo,":",xup,",",ylo,":",yup,"]"))
          sink(type="message")
          stop("Dimensions in matrix division are not equal")
        }#}}}
        #Calculate Deblended Value and assign {{{
        fmap<-wsfam/image.env$wfa[xlo:xup,ylo:yup]
        #}}}
        #Check for production of NA/NaN/Infs {{{
        if (length(which((is.na(fmap)))) != 0 ) {
          #warning("NA values produced in map. They will be zero-ed.")
          fmap[which((is.na(fmap)))]<-1.0
        }#}}}
        #}}}
      } else {
      #Aperture stamps may have errors {{{
        #If sums are equal to zero, short-cut by assigning deblended map to zero {{{
        if ((sum(wsfam)==0)|(sum(image.env$wfa[xlo:xup,ylo:yup])==0)) { fmap<-array(0, dim=dim(wsfam)) }
        #}}}
        #If sums are -ve, there may be an error {{{
        #Details {{{
        #The arrays should be +ve ± noise, or 0. Produces warning, and set to zero. }}}
        if ((sum(wsfam)<0)|(sum(image.env$wfa[xlo:xup,ylo:yup])<0)) {
          warning("The sum value of the wsfa and wfa matricies are negative
          - check for an offest error in FITS image or calculations")
          fmap<-array(0, dim=dim(wsfam))
        }#}}}
        #}}}
      }
      fmap
      #}}}
  }#}}}
  if (cutup) { registerDoParallel(cores=num.cores) } 
  #}}}

  message('=========END=======Make_Deblended_Weightmap========END=================\n')

  #Return {{{
  if (!is.null(env)) { detach(env) }
  #Return Deblended Stamps
  return=deblend.weight.map
  #}}}
}

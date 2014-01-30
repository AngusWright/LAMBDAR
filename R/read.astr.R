read.astr <-
function(fitsname){
  # Proceedure Reads the header of <fitsname> and
  # Returns the astrometry structure

  #Initialise List {{{
  astr_struc<-NULL
  #}}}
  #Check for Typical Keywords of known Lengths {{{
  astr_struc<-c(astr_struc, list(BITPIX  =read.fitskey("BITPIX",fitsname)))
  astr_struc<-c(astr_struc, list(NAXIS   =as.numeric(c(read.fitskey(c("NAXIS1","NAXIS2"),fitsname)))))
  astr_struc<-c(astr_struc, list(PC      =as.numeric(c(read.fitskey(c("PC1_1","PC1_2"),fitsname)))))
  astr_struc<-c(astr_struc, list(CDELT   =as.numeric(c(read.fitskey(c("CDELT1","CDELT2"),fitsname)))))
  astr_struc<-c(astr_struc, list(CRPIX   =as.numeric(c(read.fitskey(c("CRPIX1","CRPIX2"),fitsname)))))
  astr_struc<-c(astr_struc, list(CRVAL   =as.numeric(c(read.fitskey(c("CRVAL1","CRVAL2"),fitsname)))))
  astr_struc<-c(astr_struc, list(CD=rbind(as.numeric(c(read.fitskey(c("CD1_1","CD1_2"),fitsname))),
                                          as.numeric(c(read.fitskey(c("CD2_1","CD2_2"),fitsname))))))
  astr_struc<-c(astr_struc, list(CROTA   =as.numeric(c(read.fitskey(c("CROTA1","CROTA2"),fitsname)))))
  astr_struc<-c(astr_struc, list(LONGPOLE=as.numeric(c(read.fitskey("LONGPOLE",fitsname)))))
  astr_struc<-c(astr_struc, list(LATPOLE=as.numeric(c(read.fitskey("LATPOLE",fitsname)))))
  #}}}
  # Check for PV2 Keywords (can be up to 21) {{{
  if ((is.finite(read.fitskey("PV2_1",fitsname)))) {
    PV2<-foreach (i=1:21, .combine="cbind") %do% {
      as.numeric(read.fitskey(paste("PV2_",i,sep=""),fitsname))
    }
    PV2<-PV2[which(is.finite(PV2))]
    astr_struc<-c(astr_struc, list(PV2=PV2))
  }#}}}
  #Convert to standard Format -> use CD matrix {{{
  if (!is.finite(astr_struc$CD[1,1])) {
  astr_struc$CD<-rbind(as.numeric(c(astr_struc$CDELT[1],astr_struc$CROTA[1])),
                       as.numeric(c(astr_struc$CROTA[2],astr_struc$CDELT[2])))
  }#}}}
  #Return {{{
  return=astr_struc
  #}}}
}

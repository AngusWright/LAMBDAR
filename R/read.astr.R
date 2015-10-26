read.astr <-
function(fitsname, hdu=1){
#Proceedure Reads the header of <fitsname> and
#Returns the astrometry structure

  #Initialise List {{{
  astr_struc<-NULL
  #}}}
  oldwarn<-options(warn=-1)
  #Check for Typical Keywords of known Lengths {{{
  astr_struc<-c(astr_struc, list(BITPIX  =read.fitskey("BITPIX",fitsname,hdu=hdu)))
  astr_struc<-c(astr_struc, list(CTYPE   =as.numeric(c(read.fitskey(c("CTYPE1","CTYPE2","CTYPE3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(NAXIS   =as.numeric(c(read.fitskey(c("NAXIS1","NAXIS2","NAXIS3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(PC      =as.numeric(c(read.fitskey(c("PC1_1","PC1_2","PC1_3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(CDELT   =as.numeric(c(read.fitskey(c("CDELT1","CDELT2","CDELT3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(CRPIX   =as.numeric(c(read.fitskey(c("CRPIX1","CRPIX2","CRPIX3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(CRVAL   =as.numeric(c(read.fitskey(c("CRVAL1","CRVAL2","CRVAL3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(CD=rbind(as.numeric(c(read.fitskey(c("CD1_1","CD1_2","CD_1_3"),fitsname,hdu=hdu))),
                                          as.numeric(c(read.fitskey(c("CD2_1","CD2_2","CD_2_3"),fitsname,hdu=hdu))),
                                          as.numeric(c(read.fitskey(c("CD3_1","CD3_2","CD_3_3"),fitsname,hdu=hdu))))))
  astr_struc<-c(astr_struc, list(CROTA   =as.numeric(c(read.fitskey(c("CROTA1","CROTA2","CROTA3"),fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(LONGPOLE=as.numeric(c(read.fitskey("LONGPOLE",fitsname,hdu=hdu)))))
  astr_struc<-c(astr_struc, list(LATPOLE=as.numeric(c(read.fitskey("LATPOLE",fitsname,hdu=hdu)))))
  #}}}
  # Check for PV2 Keywords (can be up to 21) {{{
  if ((is.finite(read.fitskey("PV2_1",fitsname,hdu=hdu)))) {
    PV2<-foreach (i=1:21, .combine="cbind") %do% {
      as.numeric(read.fitskey(paste("PV2_",i,sep=""),fitsname,hdu=hdu))
    }
    PV2<-PV2[which(is.finite(PV2))]
    astr_struc<-c(astr_struc, list(PV2=PV2))
  }#}}}
  #Convert to standard Format -> use CD matrix {{{
  if (!is.finite(astr_struc$CD[1,1])) {
  astr_struc$CD<-rbind(as.numeric(c(astr_struc$CDELT[1],astr_struc$CROTA[1],astr_struc$CROTA[1])),
                       as.numeric(c(astr_struc$CROTA[2],astr_struc$CDELT[2],astr_struc$CROTA[2])),
                       as.numeric(c(astr_struc$CROTA[3],astr_struc$CROTA[3],astr_struc$CDELT[3])))
  }#}}}
  options(oldwarn)
  #Return {{{
  return(astr_struc)
  #}}}
}

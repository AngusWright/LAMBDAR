writesfatableout <-
function(env=NULL,filename) {
  #Write out a table containing the flux results for each galaxy
  # Load Parameter Space
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  #on.exit(detach('env'))

  sz=length(ra_g)

  #Get distance from image edge {{{
  #X-axis {{{
  #Distance from low side in pixels {{{
  dx_p_lo<-x_p
  #}}}
  #Distance from high side in pixels {{{
  dx_p_hi<-astr_struc$NAXIS[1]-x_p
  #}}}
  #Which apertures are closer to the high side {{{
  index_x<-which(dx_p_hi < dx_p_lo)
  #}}}
  #Distances {{{
  dx_p<-dx_p_lo
  dx_p[index_x]<-dx_p_hi[index_x]
  #}}}
  #}}}
  #Y-axis {{{
  #Distance from low side in pixels {{{
  dy_p_lo<-y_p
  #}}}
  #Distance from high side in pixels {{{
  dy_p_hi<-astr_struc$NAXIS[2]-y_p
  #}}}
  #Which apertures are closer to the high side {{{
  index_y<-which(dy_p_hi < dy_p_lo)
  #}}}
  #Distances {{{
  dy_p<-dy_p_lo
  dy_p[index_y]<-dy_p_hi[index_y]
  #}}}
  #}}}
  #}}}

  if (verbose|diagnostic) {
    newtable<-data.frame(CATA_INDEX = id_g)
    newtable[["ALPHA_J2000"]] = ra_g
    newtable[["DELTA_J2000"]] = dec_g
    newtable[["X_IMAGE"]] = x_p
    newtable[["Y_IMAGE"]] = y_p
    newtable[["NX_PIX2EDGE"]] = dx_p
    newtable[["NY_PIX2EDGE"]] = dy_p
    newtable[["THETA_J2000"]] = theta_g
    newtable[["MAJOR_ARCSEC"]] = a_g
    newtable[["MINOR_ARCSEC"]] = b_g
    newtable[["SumSA"]] = ssa
    newtable[["SumSFA"]] = ssfa
    newtable[["SumSFAsq"]] = ssfa2
    newtable[["SumSFAxData"]] = ssfad
    newtable[["SumSFAsqxErrorsq"]] = ssfa2e2
    newtable[["SkyLocal"]] = skylocal
    newtable[["SkyFlux_units"]] = skyflux
    newtable[["SkyError_units"]] = skyerr
    newtable[["SkyFlux_Jy"]] = skyflux*10^(-(magZP-8.9)/2.5)
    newtable[["SkyError_Jy"]] = (skyerr)*(10^(-(magZP-8.9)/2.5))
    newtable[["SkyRMS_units"]] = skyrms
    newtable[["SkyRMS_Jy"]] = skyrms*(10^(-(magZP-8.9)/2.5))
    newtable[["SkyRMSpval"]] = skypval
    newtable[["DetecThres_5sig"]] = detecthres
    newtable[["DetecThresMag_5sig"]] = detecthres.mag
    newtable[["SFAFlux_units"]] = sfaflux
    newtable[["SFAErr_units"]] = sfaerr
    newtable[["SFAFlux_Jy"]] = sfaflux*10^(-(magZP-8.9)/2.5)
    newtable[["SFAErr_Jy"]] = (sfaerr)*(10^(-(magZP-8.9)/2.5))
    newtable[["SumDFA"]] = sdfa
    newtable[["SumDFAsq"]] = sdfa2
    newtable[["SumDFAxData"]] = sdfad
    newtable[["SumDFAsqxErrorsq"]] = sdfa2e2
    newtable[["DFAFlux_units"]] = dfaflux
    newtable[["DFAErr_units"]] = dfaerr
    if (iterateFluxes) { for (i in 1:length(fluxiters[1,])) { newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[,i]*10^(-(magZP-8.9)/2.5) } }
    newtable[["DFAFlux_Jy"]] = dfaflux*10^(-(magZP-8.9)/2.5)
    newtable[["DFAErr_Jy"]] = (dfaerr)*(10^(-(magZP-8.9)/2.5))
    newtable[["ABMagnitude"]] = mags
    newtable[["ABMagErr"]] = (2.5/log(10))*(dfaerr/dfaflux)
    newtable[["ApCorr"]] = 1+(abs(spsf-ssfap)/spsf)
    newtable[["PixelFlux"]] = pixflux
    write.csv(newtable,file=filename, na="-", row.names=FALSE)
  } else {
    newtable=data.frame(CATA_INDEX=id_g)
    newtable[["ALPHA_J2000"]] = ra_g
    newtable[["DELTA_J2000"]] = dec_g
    newtable[["X_IMAGE"]] = x_p
    newtable[["Y_IMAGE"]] = y_p
    newtable[["NX_PIX2EDGE"]] = dx_p
    newtable[["NY_PIX2EDGE"]] = dy_p
    newtable[["SkyLocal"]] = skylocal
    newtable[["SkyFlux_units"]] = skyflux
    newtable[["SkyError_units"]] = skyerr
    newtable[["SkyFlux_Jy"]] = skyflux*(10^(-(magZP-8.9)/2.5))
    newtable[["SkyError_Jy"]] = skyerr*(10^(-(magZP-8.9)/2.5))
    newtable[["SkyRMS_units"]] = skyrms
    newtable[["SkyRMS_Jy"]] = skyrms*(10^(-(magZP-8.9)/2.5))
    newtable[["SkyRMSpval"]] = skypval
    newtable[["DetecThres_5sig"]] = detecthres
    newtable[["DetecThresMag_5sig"]] = detecthres.mag
    newtable[["SumSFA"]] = ssfa
    newtable[["SumDFA"]] = sdfa
    newtable[["SFAFlux_units"]] = sfaflux
    newtable[["SFAErr_units"]] = sfaerr
    newtable[["SFAFlux_Jy"]] = sfaflux*(10^(-(magZP-8.9)/2.5))
    newtable[["SFAErr_Jy"]] = sfaerr*(10^(-(magZP-8.9)/2.5))
    newtable[["DFAFlux_units"]] = dfaflux
    newtable[["DFAErr_units"]] = dfaerr
    if (iterateFluxes) { for (i in 1:dim(fluxiters)[2]) { newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[,i]*10^(-(magZP-8.9)/2.5) } }
    newtable[["DFAFlux_Jy"]] = dfaflux*(10^(-(magZP-8.9)/2.5))
    newtable[["DFAErr_Jy"]] = dfaerr*(10^(-(magZP-8.9)/2.5))
    newtable[["ABMagnitude"]] = mags
    newtable[["ABMagErr"]] = (2.5/log(10))*(dfaerr/dfaflux)
    newtable[["ApCorr"]] = 1+(abs(spsf-ssfap)/spsf)
    newtable[["PixelFlux"]] = pixflux
    write.csv(newtable,file=filename, na="-", row.names=FALSE)
  }

  if (!is.null(env)) { detach(env) }
  return=NULL
}

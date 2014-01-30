writesfatableout <-
function(env=NULL,filename) {
  #Write out a table containing the flux results for each galaxy
  # Load Parameter Space
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))

  sz=length(ra_g)
  newtable = t(replicate(sz, list('CATA_INDEX'=0,
                            'ALPHA_J2000'=0,
                            'DELTA_J2000'=0,
                            'X_IMAGE'=0,
                            'Y_IMAGE'=0,
                            'NX_PIX2EDGE'=0,
                            'NY_PIX2EDGE'=0,
                            'THETA_J2000'=0.,
                            'MAJOR_ARCSEC'=0,
                            'MINOR_ARCSEC'=0,
                            'SumSA'=0.,
                            'SumSFA'=0.,
                            'SumSFAsq'=0.,
                            'SumSFAxData'=0.,
                            'SumSFAsqxErrorsq'=0.,
                            'SkyLocal'=0.,
                            'SkyFlux'=0.,
                            'SkyError'=0.,
                            'SkyRMS'=0.,
                            'SkyRMSpval'=0.,
                            'DetecThres_5sig'=0.,
                            'DetecThresMag_5sig'=0.,
                            'SFAFlux'=0.,
                            'SFAErr'=0.,
                            'SumDFA'=0.,
                            'SumDFAsq'=0.,
                            'SumDFAxData'=0.,
                            'SumDFAsqxErrorsq'=0.,
                            'DFAFlux'=0.,
                            'DFAErr'=0.,
                            'Magnitude'=0.,
                            'PixelFlux'=0.)))

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

  newtable[,"CATA_INDEX"] = as.integer(id_g)
  newtable[,"ALPHA_J2000"] = ra_g
  newtable[,"DELTA_J2000"] = dec_g
  newtable[,"X_IMAGE"] = x_p
  newtable[,"Y_IMAGE"] = y_p
  newtable[,"NX_PIX2EDGE"] = dx_p
  newtable[,"NY_PIX2EDGE"] = dy_p
  newtable[,"THETA_J2000"] = theta_g
  newtable[,"MAJOR_ARCSEC"] = a_g
  newtable[,"MINOR_ARCSEC"] = b_g
  newtable[,"SumSA"] = ssa
  newtable[,"SumSFA"] = ssfa
  newtable[,"SumSFAsq"] = ssfa2
  newtable[,"SumSFAxData"] = ssfad
  newtable[,"SumSFAsqxErrorsq"] = ssfa2e2
  newtable[,"SkyLocal"] = skylocal
  newtable[,"SkyFlux"] = skyflux
  newtable[,"SkyError"] = skyerr
  newtable[,"SkyRMS"] = skyrms
  newtable[,"SkyRMSpval"] = skypval
  newtable[,"DetecThres_5sig"] = detecthres
  newtable[,"DetecThresMag_5sig"] = detecthres.mag
  newtable[,"SFAFlux"] = sfaflux
  newtable[,"SFAErr"] = sfaerr
  newtable[,"SumDFA"] = sdfa
  newtable[,"SumDFAsq"] = sdfa2
  newtable[,"SumDFAxData"] = sdfad
  newtable[,"SumDFAsqxErrorsq"] = sdfa2e2
  newtable[,"DFAFlux"] = dfaflux
  newtable[,"Magnitude"] = mags
  newtable[,"DFAErr"] = dfaerr
  newtable[,"PixelFlux"] = pixflux
  write.table(newtable,file=filename,sep=",", col.names=TRUE, na="-", dec=".", row.names=FALSE)
  detach(env)
  return=NULL
}

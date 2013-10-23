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
                            'ALPHA_J2000_R'=0,
                            'DELTA_J2000_R'=0,
                            'THETA_J2000_R'=0.,
                            'MAJOR_ARCSEC'=0,
                            'MINOR_ARCSEC'=0,
                            'SumSA'=0.,
                            'SumSFA'=0.,
                            'SumSFAsq'=0.,
                            'SumSFAxData'=0.,
                            'SumSFAsqxErrorsq'=0.,
                            'SkyFlux'=0.,
                            'SkyError'=0.,
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

  newtable[,"CATA_INDEX"] = as.integer(id_g)
  newtable[,"ALPHA_J2000_R"] = ra_g
  newtable[,"DELTA_J2000_R"] = dec_g
  newtable[,"THETA_J2000_R"] = theta_g
  newtable[,"MAJOR_ARCSEC"] = a_g
  newtable[,"MINOR_ARCSEC"] = b_g
  newtable[,"SumSA"] = ssa
  newtable[,"SumSFA"] = ssfa
  newtable[,"SumSFAsq"] = ssfa2
  newtable[,"SumSFAxData"] = ssfad
  newtable[,"SumSFAsqxErrorsq"] = ssfa2e2
  newtable[,"SkyFlux"] = skyflux
  newtable[,"SkyError"] = skyerr
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
}

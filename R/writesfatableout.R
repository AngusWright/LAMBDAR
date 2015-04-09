writesfatableout <-
function(env=NULL,filename) {
  #Write out a table containing the flux results for each galaxy
  # Load Parameter Space
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  #on.exit(detach('env'))

  sz=length(ra_g)
  janskyConv<-10^(-(magZP-8.9)/2.5)

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
    newtable[["X_IMAGE"]] = x_g
    newtable[["Y_IMAGE"]] = y_g
    newtable[["NX_PIX2EDGE"]] = dx_p
    newtable[["NY_PIX2EDGE"]] = dy_p
    newtable[["THETA_J2000"]] = theta_g
    newtable[["MAJOR_ARCSEC"]] = a_g
    newtable[["MINOR_ARCSEC"]] = b_g
    newtable[["SumSA"]] = ssa
    newtable[["SumSFA"]] = ssfa
    newtable[["SumSFAsq"]] = ssfa2
    newtable[["SumPSF"]] = spsf
    newtable[["SumSFAxPSF"]] = ssfap
    newtable[["SumSFAxData"]] = ssfad
    newtable[["SumSFAsqxErrorsq"]] = ssfa2e2
    newtable[["SkyLocal"]] = skylocal
    newtable[["SkyFlux_units"]] = skyflux
    newtable[["SkyError_units"]] = skyerr
    newtable[["SkyFlux_Jy"]] = skyflux*janskyConv
    newtable[["SkyError_Jy"]] = (skyerr)*janskyConv
    newtable[["SkyRMS_units"]] = skyrms
    newtable[["SkyRMS_Jy"]] = skyrms*janskyConv
    newtable[["SkyRMSpval"]] = skypval
    newtable[["DetecThres_5sig"]] = detecthres
    newtable[["DetecThresMag_5sig"]] = detecthres.mag
    newtable[["SFAFlux_units"]] = sfaflux
    newtable[["SFAErr_units"]] = sfaerr
    newtable[["SFAFlux_Jy"]] = sfaflux*janskyConv
    newtable[["SFAErr_Jy"]] = (sfaerr)*janskyConv
    newtable[["SumDFA"]] = sdfa
    newtable[["SumDFAsq"]] = sdfa2
    newtable[["SumDFAxData"]] = sdfad
    newtable[["SumDFAxErrorsq"]] = sdfae2
    newtable[["SumDFAsqxErrorsq"]] = sdfa2e2
    newtable[["DFAFlux_units"]] = dfaflux
    newtable[["DFAErr_units"]] = dfaerr
    if (iterateFluxes &  Magnitudes) {
      for (i in 1:length(fluxiters[1,])) {
        newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[,i]*janskyConv
        newtable[[paste("DFAErr_Iter",i,"_Jy",sep="")]] = erriters[,i]*janskyConv
        newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i]
      }
    }
    if (iterateFluxes & !Magnitudes) {
      for (i in 1:length(fluxiters[1,])) {
        newtable[[paste("DFAFlux_Iter",i,"_units",sep="")]] = fluxiters[,i]
        newtable[[paste("DFAErr_Iter",i,"_units",sep="")]] = erriters[,i]
        newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i]
      }
    }
    newtable[["DFAFlux_Jy"]] = dfaflux*janskyConv
    newtable[["DFAErr_Jy"]] = (dfaerr)*janskyConv
    newtable[["ABMagnitude"]] = mags
    newtable[["ABMagErr"]] = (2.5/log(10))*(dfaerr/dfaflux)
    newtable[["ApCorr"]] = ApCorr
    if (RanCor) {
      newtable[["RandomsMeanMean_Units"]] = randoms$randMean.mean
      newtable[["RandomsMeanSD_Units"]] = randoms$randMean.SD
      newtable[["RandomsMedMean_Units"]] = randoms$randMed.mean
      newtable[["RandomsMedSD_Units"]] = randoms$randMed.SD
      if (Magnitudes) {
        newtable[["RandomsMeanMean_Jy"]] = randoms$randMean.mean*janskyConv
        newtable[["RandomsMeanSD_Jy"]] = randoms$randMean.SD*janskyConv
        newtable[["RandomsMedMean_Jy"]] = randoms$randMed.mean*janskyConv
        newtable[["RandomsMedSD_Jy"]] = randoms$randMed.SD*janskyConv
      }
      #newtable[["Randoms_nNoMask"]] = randoms$nNoMask
      #newtable[["Randoms_meanMaskFrac"]] = randoms$meanMaskFrac
    }
    if (Magnitudes) {
      newtable[["PixelFlux_Jy"]] = pixflux*janskyConv
    } else {
      newtable[["PixelFlux_Units"]] = pixflux
    }
    newtable[["QuarteredPhot1_Units"]] = qssfad[,1]
    newtable[["QuarteredPhot2_Units"]] = qssfad[,2]
    newtable[["QuarteredPhot3_Units"]] = qssfad[,3]
    newtable[["QuarteredPhot4_Units"]] = qssfad[,4]
    newtable[["QuarteredPhot1_Deblended_Units"]] = qsdfad[,1]
    newtable[["QuarteredPhot2_Deblended_Units"]] = qsdfad[,2]
    newtable[["QuarteredPhot3_Deblended_Units"]] = qsdfad[,3]
    newtable[["QuarteredPhot4_Deblended_Units"]] = qsdfad[,4]
    #newtable[["PhotometryWarning"]] = ifelse(any(qssfad
    write.csv(newtable,file=filename, na="-", row.names=FALSE)
  } else {
    newtable=data.frame(CATA_INDEX=id_g)
    newtable[["ALPHA_J2000"]] = ra_g
    newtable[["DELTA_J2000"]] = dec_g
    newtable[["X_IMAGE"]] = x_g
    newtable[["Y_IMAGE"]] = y_g
    newtable[["NX_PIX2EDGE"]] = dx_p
    newtable[["NY_PIX2EDGE"]] = dy_p
    newtable[["SkyLocal"]] = skylocal
    newtable[["SkyFlux_units"]] = skyflux
    newtable[["SkyError_units"]] = skyerr
    newtable[["SkyFlux_Jy"]] = skyflux*janskyConv
    newtable[["SkyError_Jy"]] = skyerr*janskyConv
    newtable[["SkyRMS_units"]] = skyrms
    newtable[["SkyRMS_Jy"]] = skyrms*janskyConv
    newtable[["SkyRMSpval"]] = skypval
    newtable[["DetecThres_5sig"]] = detecthres
    newtable[["DetecThresMag_5sig"]] = detecthres.mag
    newtable[["SumSFA"]] = ssfa
    newtable[["SumDFA"]] = sdfa
    newtable[["SFAFlux_units"]] = sfaflux
    newtable[["SFAErr_units"]] = sfaerr
    newtable[["SFAFlux_Jy"]] = sfaflux*janskyConv
    newtable[["SFAErr_Jy"]] = sfaerr*janskyConv
    newtable[["DFAFlux_units"]] = dfaflux
    newtable[["DFAErr_units"]] = dfaerr
    if (iterateFluxes &  Magnitudes) {
      for (i in 1:length(fluxiters[1,])) {
        newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[,i]*janskyConv
        newtable[[paste("DFAErr_Iter",i,"_Jy",sep="")]] = erriters[,i]*janskyConv
        newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i]
      }
    }
    if (iterateFluxes & !Magnitudes) {
      for (i in 1:length(fluxiters[1,])) {
        newtable[[paste("DFAFlux_Iter",i,"_units",sep="")]] = fluxiters[,i]
        newtable[[paste("DFAErr_Iter",i,"_units",sep="")]] = erriters[,i]
        newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i]
      }
    }
    newtable[["DFAFlux_Jy"]] = dfaflux*janskyConv
    newtable[["DFAErr_Jy"]] = dfaerr*janskyConv
    newtable[["ABMagnitude"]] = mags
    newtable[["ABMagErr"]] = (2.5/log(10))*(dfaerr/dfaflux)
    newtable[["ApCorr"]] = ApCorr
    if (RanCor) {
      newtable[["RandomsMeanMean_Units"]] = randoms$randMean.mean
      newtable[["RandomsMeanSD_Units"]] = randoms$randMean.SD
      newtable[["RandomsMedMean_Units"]] = randoms$randMed.mean
      newtable[["RandomsMedSD_Units"]] = randoms$randMed.SD
      if (Magnitudes) {
        newtable[["RandomsMeanMean_Jy"]] = randoms$randMean.mean*janskyConv
        newtable[["RandomsMeanSD_Jy"]] = randoms$randMean.SD*janskyConv
        newtable[["RandomsMedMean_Jy"]] = randoms$randMed.mean*janskyConv
        newtable[["RandomsMedSD_Jy"]] = randoms$randMed.SD*janskyConv
      }
      #newtable[["Randoms_nNoMask"]] = randoms$nNoMask
      #newtable[["Randoms_meanMaskFrac"]] = randoms$meanMaskFrac
    }
    if (Magnitudes) {
      newtable[["PixelFlux_Jy"]] = pixflux*janskyConv
    } else {
      newtable[["PixelFlux_Units"]] = pixflux
    }
    newtable[["QuarteredPhot1_Units"]] = qssfad[,1]
    newtable[["QuarteredPhot2_Units"]] = qssfad[,2]
    newtable[["QuarteredPhot3_Units"]] = qssfad[,3]
    newtable[["QuarteredPhot4_Units"]] = qssfad[,4]
    newtable[["QuarteredPhot1_Deblended_Units"]] = qsdfad[,1]
    newtable[["QuarteredPhot2_Deblended_Units"]] = qsdfad[,2]
    newtable[["QuarteredPhot3_Deblended_Units"]] = qsdfad[,3]
    newtable[["QuarteredPhot4_Deblended_Units"]] = qsdfad[,4]
    #newtable[["PhotometryWarning"]] = ifelse(any(qssfad
    write.csv(newtable,file=filename, na="-", row.names=FALSE)
  }

  if (!is.null(env)) { detach(env) }
  return=NULL
}

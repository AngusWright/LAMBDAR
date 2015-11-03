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
  notes<-sub(".csv","_notes.txt",filename)
  sink(type="output",file=notes)
  on.exit(sink(type='output'),add=TRUE)
  today <- Sys.Date()

  cat(paste0("# Notes file for LAMBDAR (v ",packageVersion('LAMBDAR'),") Catalogue created ", format(today, format="%B %d %Y"),"\n"))
  cat(paste0("# Below lists all parameters in the accompaning catalogue, and what they represent.\n"))
  if (!PSFWeighted) {
    cat(paste0("# NB: This run of LAMBDAR did not request PSF Weighted Apertures. As such, each aperture was passed through a binary-filter which\n",
        "#     converted the contiguous apertures into TopHat apertures. The filter point was requested at ", paste(apLimit),", which means that\n",
        "#     pixels contained within the ", paste(apLimit*100),"% contour was set to 1, and the rest were set to 0.\n",
        "# IMPORTANT: Because of this, all references below to the 'Post-Convolution Aperture' (i.e. SFA & DFA) are also after the aperture was binary-filtered.\n"))
  }
  if (Magnitudes) {
    cat(paste0("# As Magnitudes were requested, all values below have been provided converted to Jansky from image units.\n",
        "# The image-to-jansky conversion factor is ",janskyConv,"\n"))
  }

  if (verbose|diagnostic) {
####CATALOGUE PARAMETER
    newtable<-data.frame(CATA_INDEX = id_g)
    colnames(newtable)<-catalab; cat(paste0(catalab," #Catalogue ID. Duplicates are prepended with 'DuplicatedID_'\n"))
    newtable[[ralab]] = ra_g; cat(paste0(ralab," #Right Ascention. Same as provided in the input catalogue\n"))
    newtable[[declab]] = dec_g; cat(paste0(declab," #Declination. Same as provided in the input catalogue\n"))
    newtable[[thetalab]] = theta_g; cat(paste0(thetalab," #Aperture Orientation Angle. Same as provided in the input catalogue\n"))
    newtable[[semimajlab]] = a_g; cat(paste0(semimajlab," #Aperture SemiMajor Axis length, in arcseconds. Same as provided in the input catalogue\n"))
    newtable[[semiminlab]] = b_g; cat(paste0(semiminlab," #Aperture SemiMinor Axis length, in arcseconds. Same as provided in the input catalogue\n"))
    if (exists("contam")) { newtable[[contamlab]] = contams; cat(paste0(contamlab," #Contaminant Flag. Same as provided in the input catalogue\n")) }
#   IMAGE PARAMETERS
    newtable[["X_IMAGE"]] = x_g; cat(paste0("X_IMAGE"," #Aperture Location in Image Coordinates; X-axis\n"))
    newtable[["Y_IMAGE"]] = y_g; cat(paste0("Y_IMAGE"," #Aperture Location in Image Coordinates; Y-axis\n"))
    newtable[["NX_PIX2EDGE"]] = dx_p; cat(paste0("NX_PIX2EDGE"," #Number of Pixels to the Nearest Edge of the Image; X-Axis. Useful for avoiding selection effects at image edges\n"))
    newtable[["NY_PIX2EDGE"]] = dy_p; cat(paste0("NY_PIX2EDGE"," #Number of Pixels to the Nearest Edge of the Image; Y-Axis. Useful for avoiding selection effects at image edges\n"))
    newtable[["StampDiameter_pix"]] = stamplen; cat(paste0("StampDiameter_pix"," #Diameter of the Aperture Stamp, in Pixels.\n"))
#   APERTURE SUMS & INTERMEDIATE SUMS
    newtable[["SumPSF"]] = spsf; cat(paste0("SumPSF"," #Integral of the PSF. Should be the same for every object\n"))
    newtable[["SumSA"]] = ssa; cat(paste0("SumSA"," #Integral of the Catalogue Aperture (i.e. the 'Prior') in pixels\n"))
    newtable[["SumSFA"]] = ssfa; cat(paste0("SumSFA"," #Integral of the Post-Convolution Aperture.\n"))
    newtable[["SumDFA"]] = sdfa; cat(paste0("SumDFA"," #Integral of the Deblended Post-Convolution Aperture.\n"))
    newtable[["SumSFAsq"]] = ssfa2; cat(paste0("SumSFAsq"," #Integral of the Post-Convolution Aperture squared.\n"))
    newtable[["SumDFAsq"]] = sdfa2; cat(paste0("SumDFAsq"," #Integral of the Deblended Post-Convolution Aperture squared.\n"))
    newtable[["SumSFAxPSF"]] = ssfap; cat(paste0("SumSFAxPSF"," #Integral of the Post-Convolution Aperture times the PSF. Used in calculating the Minimum Aperture Correction.\n"))
    newtable[["SumSFAxData"]] = ssfad; cat(paste0("SumSFAxData"," #Integral of the Post-Convolution Aperture times the Data Image. The simplest flux LAMBDAR provides.\n"))
    newtable[["SumDFAxData"]] = sdfad; cat(paste0("SumDFAxData"," #Integral of the Deblended Post-Convolution Aperture times the Data Image. The raw deblended flux LAMBDAR.\n"))
    newtable[["SumSFAxError"]] = ssfae; cat(paste0("SumSFAxError"," #Integral of the Post-Convolution Aperture times the Error Image.\n"))
    newtable[["SumDFAxError"]] = sdfae; cat(paste0("SumDFAxError"," #Integral of the Deblended Post-Convolution Aperture times the Error Image\n"))
    newtable[["SumSFAxErrorsq"]] = ssfae2; cat(paste0("SumSFAxErrorsq"," #Integral of the Post-Convolution Aperture times the Error Image squared\n"))
    newtable[["SumDFAxErrorsq"]] = sdfae2; cat(paste0("SumDFAxErrorsq"," #Integral of the Deblended Post-Convolution Aperture times the Error Image squared\n"))
    newtable[["SumSFAsqxErrorsq"]] = ssfa2e2; cat(paste0("SumSFAsqxErrorsq"," #Integral of the Post-Convolution Aperture squared times the Error Image squared\n"))
    newtable[["SumDFAsqxErrorsq"]] = sdfa2e2; cat(paste0("SumDFAsqxErrorsq"," #Integral of the Deblended Post-Convolution Aperture squared times the Error Image squared\n"))
    if (doskyest | getskyrms) {
      if (Magnitudes) {
#       SKY PARAMETERS: Jansky
        newtable[["SkyLocal_Jy"]] = skylocal*janskyConv; cat(paste0("SkyLocal_Jy"," #The measured Median Sky Level, in Jansky.\n"))
        newtable[["SkyLocal_Mean_Jy"]] = skylocal.mean*janskyConv; cat(paste0("SkyLocal_Mean_Jy"," #The measured Mean Sky Level, in Jansky.\n"))
        newtable[["SkyFlux_Jy"]] = skyflux*janskyConv; cat(paste0("SkyFlux_Jy"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Uses the measured Median Sky Level.\n"))
        newtable[["SkyFlux_Mean_Jy"]] = skyflux.mean*janskyConv; cat(paste0("SkyFlux_Mean_Jy"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Uses the measured Mean Sky Level.\n"))
        newtable[["SkyError_Jy"]] = skyerr*janskyConv; cat(paste0("SkyError_Jy"," #Error on the Median Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Determined by measuring the Error on the Median Sky Level.\n"))
        newtable[["SkyError_Mean_Jy"]] = skyerr.mean*janskyConv; cat(paste0("SkyError_Mean_Jy"," #Error on the Mean Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Determined by measuring the Error on the Mean Sky Level.\n"))
        newtable[["SkyRMS_Jy"]] = skyrms*janskyConv; cat(paste0("SkyRMS_Jy"," #RMS of the Sky Pixels calculated during the Median Sky Level calculation, in Jansky.\n"))
        newtable[["SkyRMS_Mean_Jy"]] = skyrms.mean*janskyConv; cat(paste0("SkyRMS_Mean_Jy"," #RMS of the Sky Pixels calculated during the Mean Sky Level calculation, in Jansky.\n"))
      } else {
#       SKY PARAMETERS: units
        newtable[["SkyLocal_units"]] = skylocal; cat(paste0("SkyLocal_units"," #The measured Median Sky Level, in image units.\n"))
        newtable[["SkyLocal_Mean_units"]] = skylocal.mean; cat(paste0("SkyLocal_Mean_units"," #The measured Mean Sky Level, in image units.\n"))
        newtable[["SkyFlux_units"]] = skyflux; cat(paste0("SkyFlux_units"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Uses the measured Median Sky Level.\n"))
        newtable[["SkyFlux_Mean_units"]] = skyflux.mean; cat(paste0("SkyFlux_Mean_units"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Uses the measured Mean Sky Level.\n"))
        newtable[["SkyError_units"]] = skyerr; cat(paste0("SkyError_units"," #Error on the Median Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Determined by measuring the Error on the Median Sky Level.\n"))
        newtable[["SkyError_Mean_units"]] = skyerr.mean; cat(paste0("SkyError_Mean_units"," #Error on the Mean Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Determined by measuring the Error on the Mean Sky Level.\n"))
        newtable[["SkyRMS_units"]] = skyrms; cat(paste0("SkyRMS_units"," #RMS of the Sky Pixels calculated during the Median Sky Level calculation, in image units.\n"))
        newtable[["SkyRMS_Mean_units"]] = skyrms.mean; cat(paste0("SkyRMS_Mean_units"," #RMS of the Sky Pixels calculated during the Mean Sky Level calculation, in image units.\n"))
      }
#     SKY PARAMETERS: General
      newtable[["SkyRMSpval"]] = skypval; cat(paste0("SkyRMSpval"," #p-value resulting from a Pearson's Test of Normality on Sky Pixels used in Median Sky Level calculations. Useful in determining behaviour of sky pixels\n"))
      newtable[["SkyRMSpval_Mean"]] = skypval.mean; cat(paste0("SkyRMSpval_Mean"," #p-value resulting from a Pearson's Test of Normality on Sky Pixels used in Mean Sky Level calculations. Useful in determining behaviour of sky pixels\n"))
      newtable[["NearSkyBins"]] = skyNBinNear; cat(paste0("NearSkyBins"," #Number of Sky Bins equal to the Sky Estimate (within uncertainties), during Median Sky Level Calculations. Max 10.\n"))
      newtable[["NearSkyBins_Mean"]] = skyNBinNear.mean; cat(paste0("NearSkyBins_Mean"," #Number of Sky Bins equal to the Sky Estimate (within uncertainties), during Mean Sky Level Calculations. Max 10.\n"))
    }
    if (!Magnitudes) {
#   FLUXES: Units
      newtable[["SFAFlux_units"]] = sfaflux; cat(paste0("SFAFlux_units"," #Final Blended Flux, in image units.\n"))
      newtable[["SFAErr_units"]] = sfaerr; cat(paste0("SFAErr_units"," #Error on Final Blended Flux, in image units.\n"))
      newtable[["DFAFlux_units"]] = dfaflux; cat(paste0("DFAFlux_units"," #Final Deblended Flux, in image units.\n"))
      newtable[["DFAErr_units"]] = dfaerr; cat(paste0("DFAErr_units"," #Error on Final Deblended Flux, in image units.\n"))
      newtable[["DeblendErr_units"]] = deblerr; cat(paste0("DeblendErr_units"," #Absolute Flux Uncertainty caused by the Deblend of the Deblended Post-Convolution Aperture.\n"))
      if (iterateFluxes) {
        cat(paste0(paste("DFAFlux_Iter<i>_units",sep="")," #Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. NB: These values are *not* sky subtracted and/or aperture corrected. In image units.\n"))
        cat(paste0(paste("DFAErr_Iter<i>_units",sep="")," #Error on Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. In image units.\n"))
        cat(paste0(paste("SumDFA_Iter<i>",sep="")," #Sum of the Deblended Post-Convolution Model after <i> iterations of Measurement, Fluxweighting, and Deblending. In pixels\n"))
        if (length(dim(fluxiters))>1) {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_units",sep="")]] = fluxiters[,i];
            newtable[[paste("DFAErr_Iter",i,"_units",sep="")]] = erriters[,i];
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i];
          }
        } else {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_units",sep="")]] = fluxiters[i];
            newtable[[paste("DFAErr_Iter",i,"_units",sep="")]] = erriters[i];
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[i];
          }
        }
      }
    } else {
#     FLUXES: Jansky
      newtable[["SFAFlux_Jy"]] = sfaflux*janskyConv; cat(paste0("SFAFlux_Jy"," #Final Blended Flux, in Jansky.\n"))
      newtable[["SFAErr_Jy"]] = sfaerr*janskyConv; cat(paste0("SFAErr_Jy"," #Error on Final Blended Flux, in Jansky.\n"))
      newtable[["DFAFlux_Jy"]] = dfaflux*janskyConv; cat(paste0("DFAFlux_Jy"," #Final Deblended Flux, in Jansky.\n"))
      newtable[["DFAErr_Jy"]] = dfaerr*janskyConv; cat(paste0("DFAErr_Jy"," #Error on Final Deblended Flux, in Jansky.\n"))
      newtable[["DeblendErr_Jy"]] = deblerr*janskyConv; cat(paste0("DeblendErr_Jy"," #Absolute Flux Uncertainty caused by the Deblend of the Deblended Post-Convolution Aperture, in Jansky.\n"))
      if (iterateFluxes) {
        cat(paste0(paste("DFAFlux_Iter<i>_Jy",sep="")," #Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. NB: These values are *not* sky subtracted and/or aperture corrected. In Jansky.\n"))
        cat(paste0(paste("DFAErr_Iter<i>_Jy",sep="")," #Error on Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. In Jansky.\n"))
        cat(paste0(paste("SumDFA_Iter<i>",sep="")," #Sum of the Deblended Post-Convolution Model after <i> iterations of Measurement, Fluxweighting, and Deblending. In pixels\n"))
        if (length(dim(fluxiters))>1) {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[,i]*janskyConv;
            newtable[[paste("DFAErr_Iter",i,"_Jy",sep="")]] = erriters[,i]*janskyConv;
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i];
          }
        } else {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[i]*janskyConv;
            newtable[[paste("DFAErr_Iter",i,"_Jy",sep="")]] = erriters[i]*janskyConv;
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[i];
          }
        }
      }
      newtable[["ABMagDFA"]] = mags; cat(paste0("ABMagDFA"," #AB Magnitude of Final Deblended Flux\n"))
      newtable[["ABMagErrDFA"]] = (2.5/log(10))*(dfaerr/dfaflux); cat(paste0("ABMagErrDFA"," #Error on the AB Magnitude of Final Deblended Flux\n"))
    }
#   BLANKS CORRECTION
    if (BlankCor) {
      if (Magnitudes) {
        newtable[["BlanksMeanMean_Jy"]] = blanks$randMean.mean*janskyConv; cat(paste0("BlanksMeanMean_Jy"," #Mean Value of the Mean Pixel Values in <n> Blank Apertures. In Jansky.\n"))
        newtable[["BlanksMeanSD_Jy"]] = blanks$randMean.SD*janskyConv; cat(paste0("BlanksMeanSD_Jy"," #Standard Deviation of the Mean Pixel Values in <n> Blank Apertures. In Jansky.\n"))
        newtable[["BlanksMeanMAD_Jy"]] = blanks$randMean.MAD*janskyConv; cat(paste0("BlanksMeanMAD_Jy"," #Median Absolute Deviation from Median of the Mean Pixel Values in <n> Blank Apertures. In Jansky.\n"))
        newtable[["BlanksApMean"]] = blanks$randAp.mean; cat(paste0("BlanksApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApSD"]] = blanks$randAp.SD; cat(paste0("BlanksApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApMAD"]] = blanks$randAp.MAD; cat(paste0("BlanksApMAD"," #Median Absolute Deviation from Median of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
      } else {
        newtable[["BlanksMeanMean_units"]] = blanks$randMean.mean; cat(paste0("BlanksMeanMean_units"," #Mean Value of the Mean Pixel Values in <n> Blank Apertures. In image units.\n"))
        newtable[["BlanksMeanSD_units"]] = blanks$randMean.SD; cat(paste0("BlanksMeanSD_units"," #Standard Deviation of the Mean Pixel Values in <n> Blank Apertures. In image units.\n"))
        newtable[["BlanksMeanMAD_units"]] = blanks$randMean.MAD; cat(paste0("BlanksMeanMAD_Jy"," #Median Absolute Deviation from Median of the Mean Pixel Values in <n> Blank Apertures.\n"))
        newtable[["BlanksApMean"]] = blanks$randAp.mean; cat(paste0("BlanksApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApSD"]] = blanks$randAp.SD; cat(paste0("BlanksApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApMAD"]] = blanks$randAp.MAD; cat(paste0("BlanksApMAD"," #Median Absolute Deviation from Median of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
      }
      newtable[["NumberMeasuredBlanks"]] = blanks$nRand; cat(paste0("NumberMeasuredBlanks"," #Number of Blank Apertures measured per object. Varies object to object due to source masking.\n"))
    }
    if (RanCor) {
      if (Magnitudes) {
        newtable[["RandomsMeanMean_Jy"]] = randoms$randMean.mean*janskyConv; cat(paste0("RandomsMeanMean_Jy"," #Mean Value of the Mean Pixel Values in <n> Random Apertures. In Jansky.\n"))
        newtable[["RandomsMeanSD_Jy"]] = randoms$randMean.SD*janskyConv; cat(paste0("RandomsMeanSD_Jy"," #Standard Deviation of the Mean Pixel Values in <n> Random Apertures. In Jansky.\n"))
        newtable[["RandomsApMean"]] = randoms$randAp.mean; cat(paste0("RandomsApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
        newtable[["RandomsApSD"]] = randoms$randAp.SD; cat(paste0("RandomsApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
      } else {
        newtable[["RandomsMeanMean_units"]] = randoms$randMean.mean; cat(paste0("RandomsMeanMean_units"," #Mean Value of the Mean Pixel Values in <n> Random Apertures. In image units.\n"))
        newtable[["RandomsMeanSD_units"]] = randoms$randMean.SD; cat(paste0("RandomsMeanSD_units"," #Standard Deviation of the Mean Pixel Values in <n> Random Apertures. In image units.\n"))
        newtable[["RandomsApMean"]] = randoms$randAp.mean; cat(paste0("RandomsApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
        newtable[["RandomsApSD"]] = randoms$randAp.SD; cat(paste0("RandomsApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
      }
    }
    if (!nopsf) {
      newtable[["MinApCorr"]] = ApCorr; cat(paste0("MinApCorr"," #Value of the Multiplicative Minimum Aperture Correction. This correction is important, and can _only_ improve an objects flux determination.\n"))
      newtable[["MaxApCorr"]] = WtCorr; cat(paste0("MaxApCorr"," #Value of the Multiplicative Maximum Aperture Correction. This correction can cause excess flux-boosting systematically as a function of object morphology, and is provided __as a guide only__\n"))
    }
#   Misc.
    if (Magnitudes) {
      newtable[["DetecThres_5sig_Jy"]] = detecthres*janskyConv; cat(paste0("DetecThres_5sig"," #Flux contained in this aperture, assuming a 5-sigma Detection, in Jansky\n"))
      newtable[["DetecThres_5sig_Mag"]] = detecthres.mag; cat(paste0("DetecThres_5sig_Mag"," #AB Magnitude of an object in this aperture assimung a 5-sigma Detection.\n"))
    } else {
      newtable[["DetecThres_5sig_units"]] = detecthres; cat(paste0("DetecThres_5sig"," #Flux contained in this aperture, assuming a 5-sigma Detection, in image units\n"))
    }
    if (Magnitudes) {
      cat(paste0("QuarteredPhot1_Jy"," #Flux contained in Four Quadrants of the Post-Convolution aperture, split along the stamp centre in x & y directions, in Jansky. Useful in diagnosis of bad apertures/deblends.\n"))
      cat(paste0("QuarteredPhot1_Deblended_Jy"," #Flux contained in Four Quadrants of the Deblended Post-Convolution aperture, split along the stamp centre in x & y directions, in Jansky. Useful in diagnosis of bad apertures/deblends.\n"))
      if (length(dim(qssfad))>1) {
        newtable[["QuarteredPhot1_Jy"]] = qssfad[,1]*janskyConv;
        newtable[["QuarteredPhot2_Jy"]] = qssfad[,2]*janskyConv;
        newtable[["QuarteredPhot3_Jy"]] = qssfad[,3]*janskyConv;
        newtable[["QuarteredPhot4_Jy"]] = qssfad[,4]*janskyConv;
        newtable[["QuarteredPhot1_Deblended_Jy"]] = qsdfad[,1]*janskyConv;
        newtable[["QuarteredPhot2_Deblended_Jy"]] = qsdfad[,2]*janskyConv;
        newtable[["QuarteredPhot3_Deblended_Jy"]] = qsdfad[,3]*janskyConv;
        newtable[["QuarteredPhot4_Deblended_Jy"]] = qsdfad[,4]*janskyConv;
      } else {
        newtable[["QuarteredPhot1_Jy"]] = qssfad[1]*janskyConv;
        newtable[["QuarteredPhot2_Jy"]] = qssfad[2]*janskyConv;
        newtable[["QuarteredPhot3_Jy"]] = qssfad[3]*janskyConv;
        newtable[["QuarteredPhot4_Jy"]] = qssfad[4]*janskyConv;
        newtable[["QuarteredPhot1_Deblended_Jy"]] = qsdfad[1]*janskyConv;
        newtable[["QuarteredPhot2_Deblended_Jy"]] = qsdfad[2]*janskyConv;
        newtable[["QuarteredPhot3_Deblended_Jy"]] = qsdfad[3]*janskyConv;
        newtable[["QuarteredPhot4_Deblended_Jy"]] = qsdfad[4]*janskyConv;
      }
      newtable[["PixelFlux_Jy"]] = pixflux*janskyConv; cat(paste0("PixelFlux_Jy"," #Pixel Flux at the object RA/DEC, in Jansky.\n"))
    } else {
      cat(paste0("QuarteredPhot1_units"," #Flux contained in Four Quadrants of the Post-Convolution aperture, split along the stamp centre in x & y directions, in image units. Useful in diagnosis of bad apertures/deblends.\n"))
      cat(paste0("QuarteredPhot1_Deblended_units"," #Flux contained in Four Quadrants of the Deblended Post-Convolution aperture, split along the stamp centre in x & y directions, in image units. Useful in diagnosis of bad apertures/deblends.\n"))
      if (length(dim(qssfad))>1) {
        newtable[["QuarteredPhot1_units"]] = qssfad[,1];
        newtable[["QuarteredPhot2_units"]] = qssfad[,2];
        newtable[["QuarteredPhot3_units"]] = qssfad[,3];
        newtable[["QuarteredPhot4_units"]] = qssfad[,4];
        newtable[["QuarteredPhot1_Deblended_units"]] = qsdfad[,1];
        newtable[["QuarteredPhot2_Deblended_units"]] = qsdfad[,2];
        newtable[["QuarteredPhot3_Deblended_units"]] = qsdfad[,3];
        newtable[["QuarteredPhot4_Deblended_units"]] = qsdfad[,4];
      } else {
        newtable[["QuarteredPhot1_units"]] = qssfad[1];
        newtable[["QuarteredPhot2_units"]] = qssfad[2];
        newtable[["QuarteredPhot3_units"]] = qssfad[3];
        newtable[["QuarteredPhot4_units"]] = qssfad[4];
        newtable[["QuarteredPhot1_Deblended_units"]] = qsdfad[1];
        newtable[["QuarteredPhot2_Deblended_units"]] = qsdfad[2];
        newtable[["QuarteredPhot3_Deblended_units"]] = qsdfad[3];
        newtable[["QuarteredPhot4_Deblended_units"]] = qsdfad[4];
      }
      newtable[["PixelFlux_units"]] = pixflux; cat(paste0("PixelFlux_units"," #Pixel Flux at the object RA/DEC, in image units\n"))
    }
    newtable[["PhotometryWarning"]] = photWarnFlag; cat(paste0(""," #Some Rudimentary Photometry Warnings. Read Carefully and Take Note!\n"))
    cat(paste0(""," #PhotometryWarnings: Q - Quartered Photometry Warning; 70+% of source flux contained in a single quadrant\n"))
    cat(paste0("","                      S - Sky Estimate Warning; Sky estimate determined in by 3 or less bins\n"))
    cat(paste0("","                      X - Saturation Warning; Pixels in this object aperture are saturated\n"))
    cat(paste0("","                      I - Iteration Warning; This object was lost during iteration because its flux was measured to be <= 0.\n"))
    write.csv(newtable,file=filename, na="-", row.names=FALSE)
  } else {
####CATALOGUE PARAMETER
    newtable<-data.frame(CATA_INDEX = id_g)
    colnames(newtable)<-catalab; cat(paste0(catalab," #Catalogue ID. Duplicates are prepended with 'DuplicatedID_'\n"))
    newtable[[ralab]] = ra_g; cat(paste0(ralab," #Right Ascention. Same as provided in the input catalogue\n"))
    newtable[[declab]] = dec_g; cat(paste0(declab," #Declination. Same as provided in the input catalogue\n"))
    newtable[[thetalab]] = theta_g; cat(paste0(thetalab," #Aperture Orientation Angle. Same as provided in the input catalogue\n"))
    newtable[[semimajlab]] = a_g; cat(paste0(semimajlab," #Aperture SemiMajor Axis length, in arcseconds. Same as provided in the input catalogue\n"))
    newtable[[semiminlab]] = b_g; cat(paste0(semiminlab," #Aperture SemiMinor Axis length, in arcseconds. Same as provided in the input catalogue\n"))
    if (exists("contam")) { newtable[[contamlab]] = contams; cat(paste0(contamlab," #Contaminant Flag. Same as provided in the input catalogue\n")) }
#   IMAGE PARAMETERS
    newtable[["NX_PIX2EDGE"]] = dx_p; cat(paste0("NX_PIX2EDGE"," #Number of Pixels to the Nearest Edge of the Image; X-Axis. Useful for avoiding selection effects at image edges\n"))
    newtable[["NY_PIX2EDGE"]] = dy_p; cat(paste0("NY_PIX2EDGE"," #Number of Pixels to the Nearest Edge of the Image; Y-Axis. Useful for avoiding selection effects at image edges\n"))
#   APERTURE SUMS & INTERMEDIATE SUMS
    newtable[["SumPSF"]] = spsf; cat(paste0("SumPSF"," #Integral of the PSF. Should be the same for every object\n"))
    newtable[["SumSA"]] = ssa; cat(paste0("SumSA"," #Integral of the Catalogue Aperture (i.e. the 'Prior') in pixels\n"))
    newtable[["SumSFA"]] = ssfa; cat(paste0("SumSFA"," #Integral of the Post-Convolution Aperture.\n"))
    newtable[["SumDFA"]] = sdfa; cat(paste0("SumDFA"," #Integral of the Deblended Post-Convolution Aperture.\n"))
    if (doskyest | getskyrms) {
      if (Magnitudes) {
#       SKY PARAMETERS: Jansky
        newtable[["SkyLocal_Jy"]] = skylocal*janskyConv; cat(paste0("SkyLocal_Jy"," #The measured Median Sky Level, in Jansky.\n"))
        newtable[["SkyLocal_Mean_Jy"]] = skylocal.mean*janskyConv; cat(paste0("SkyLocal_Mean_Jy"," #The measured Mean Sky Level, in Jansky.\n"))
        newtable[["SkyFlux_Jy"]] = skyflux*janskyConv; cat(paste0("SkyFlux_Jy"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Uses the measured Median Sky Level.\n"))
        newtable[["SkyFlux_Mean_Jy"]] = skyflux.mean*janskyConv; cat(paste0("SkyFlux_Mean_Jy"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Uses the measured Mean Sky Level.\n"))
        newtable[["SkyError_Jy"]] = skyerr*janskyConv; cat(paste0("SkyError_Jy"," #Error on the Median Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Determined by measuring the Error on the Median Sky Level.\n"))
        newtable[["SkyError_Mean_Jy"]] = skyerr.mean*janskyConv; cat(paste0("SkyError_Mean_Jy"," #Error on the Mean Sky Flux contained in the Deblended Post-Convolution Aperture, in Jansky. Determined by measuring the Error on the Mean Sky Level.\n"))
        newtable[["SkyRMS_Jy"]] = skyrms*janskyConv; cat(paste0("SkyRMS_Jy"," #RMS of the Sky Pixels calculated during the Median Sky Level calculation, in Jansky.\n"))
        newtable[["SkyRMS_Mean_Jy"]] = skyrms.mean*janskyConv; cat(paste0("SkyRMS_Mean_Jy"," #RMS of the Sky Pixels calculated during the Mean Sky Level calculation, in Jansky.\n"))
      } else {
#       SKY PARAMETERS: units
        newtable[["SkyLocal_units"]] = skylocal*janskyConv; cat(paste0("SkyLocal_units"," #The measured Median Sky Level, in image units.\n"))
        newtable[["SkyLocal_Mean_units"]] = skylocal.mean*janskyConv; cat(paste0("SkyLocal_Mean_units"," #The measured Mean Sky Level, in image units.\n"))
        newtable[["SkyFlux_units"]] = skyflux*janskyConv; cat(paste0("SkyFlux_units"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Uses the measured Median Sky Level.\n"))
        newtable[["SkyFlux_Mean_units"]] = skyflux.mean*janskyConv; cat(paste0("SkyFlux_Mean_units"," #Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Uses the measured Mean Sky Level.\n"))
        newtable[["SkyError_units"]] = skyerr*janskyConv; cat(paste0("SkyError_units"," #Error on the Median Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Determined by measuring the Error on the Median Sky Level.\n"))
        newtable[["SkyError_Mean_units"]] = skyerr.mean*janskyConv; cat(paste0("SkyError_Mean_units"," #Error on the Mean Sky Flux contained in the Deblended Post-Convolution Aperture, in image units. Determined by measuring the Error on the Mean Sky Level.\n"))
        newtable[["SkyRMS_units"]] = skyrms*janskyConv; cat(paste0("SkyRMS_units"," #RMS of the Sky Pixels calculated during the Median Sky Level calculation, in image units.\n"))
        newtable[["SkyRMS_Mean_units"]] = skyrms.mean*janskyConv; cat(paste0("SkyRMS_Mean_units"," #RMS of the Sky Pixels calculated during the Mean Sky Level calculation, in image units.\n"))
      }
    }
    if (!Magnitudes) {
#   FLUXES: Units
      newtable[["SFAFlux_units"]] = sfaflux; cat(paste0("SFAFlux_units"," #Final Blended Flux, in image units.\n"))
      newtable[["SFAErr_units"]] = sfaerr; cat(paste0("SFAErr_units"," #Error on Final Blended Flux, in image units.\n"))
      newtable[["DFAFlux_units"]] = dfaflux; cat(paste0("DFAFlux_units"," #Final Deblended Flux, in image units.\n"))
      newtable[["DFAErr_units"]] = dfaerr; cat(paste0("DFAErr_units"," #Error on Final Deblended Flux, in image units.\n"))
      newtable[["DeblendErr_units"]] = deblerr; cat(paste0("DeblendErr_units"," #Absolute Flux Uncertainty caused by the Deblend of the Deblended Post-Convolution Aperture.\n"))
      if (iterateFluxes) {
        cat(paste0(paste("DFAFlux_Iter<i>_units",sep="")," #Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. NB: These values are *not* sky subtracted and/or aperture corrected. In image units.\n"))
        cat(paste0(paste("DFAErr_Iter<i>_units",sep="")," #Error on Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. In image units.\n"))
        cat(paste0(paste("SumDFA_Iter<i>",sep="")," #Sum of the Deblended Post-Convolution Model after <i> iterations of Measurement, Fluxweighting, and Deblending. In pixels\n"))
        if (length(dim(fluxiters))>1) {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_units",sep="")]] = fluxiters[,i];
            newtable[[paste("DFAErr_Iter",i,"_units",sep="")]] = erriters[,i];
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i];
          }
        } else {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_units",sep="")]] = fluxiters[i];
            newtable[[paste("DFAErr_Iter",i,"_units",sep="")]] = erriters[i];
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[i];
          }
        }
      }
    } else {
#     FLUXES: Jansky
      newtable[["SFAFlux_Jy"]] = sfaflux*janskyConv; cat(paste0("SFAFlux_Jy"," #Final Blended Flux, in Jansky.\n"))
      newtable[["SFAErr_Jy"]] = sfaerr*janskyConv; cat(paste0("SFAErr_Jy"," #Error on Final Blended Flux, in Jansky.\n"))
      newtable[["DFAFlux_Jy"]] = dfaflux*janskyConv; cat(paste0("DFAFlux_Jy"," #Final Deblended Flux, in Jansky.\n"))
      newtable[["DFAErr_Jy"]] = dfaerr*janskyConv; cat(paste0("DFAErr_Jy"," #Error on Final Deblended Flux, in Jansky.\n"))
      newtable[["DeblendErr_Jy"]] = deblerr*janskyConv; cat(paste0("DeblendErr_Jy"," #Absolute Flux Uncertainty caused by the Deblend of the Deblended Post-Convolution Aperture, in Jansky.\n"))
      if (iterateFluxes) {
        cat(paste0(paste("DFAFlux_Iter<i>_Jy",sep="")," #Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. NB: These values are *not* sky subtracted and/or aperture corrected. In Jansky.\n"))
        cat(paste0(paste("DFAErr_Iter<i>_Jy",sep="")," #Error on Deblended Flux after <i> iterations of Measurement, Fluxweighting, and Deblending. In Jansky.\n"))
        cat(paste0(paste("SumDFA_Iter<i>",sep="")," #Sum of the Deblended Post-Convolution Model after <i> iterations of Measurement, Fluxweighting, and Deblending. In pixels\n"))
        if (length(dim(fluxiters))>1) {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[,i]*janskyConv;
            newtable[[paste("DFAErr_Iter",i,"_Jy",sep="")]] = erriters[,i]*janskyConv;
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[,i];
          }
        } else {
          for (i in 1:nIterations) {
            newtable[[paste("DFAFlux_Iter",i,"_Jy",sep="")]] = fluxiters[i]*janskyConv;
            newtable[[paste("DFAErr_Iter",i,"_Jy",sep="")]] = erriters[i]*janskyConv;
            newtable[[paste("SumDFA_Iter",i,sep="")]] = sdfaiters[i];
          }
        }
      }
      newtable[["ABMagDFA"]] = mags; cat(paste0("ABMagDFA"," #AB Magnitude of Final Deblended Flux\n"))
      newtable[["ABMagErrDFA"]] = (2.5/log(10))*(dfaerr/dfaflux); cat(paste0("ABMagErrDFA"," #Error on the AB Magnitude of Final Deblended Flux\n"))
    }
#   BLANKS CORRECTION
    if (BlankCor) {
      if (Magnitudes) {
        newtable[["BlanksMeanMean_Jy"]] = blanks$randMean.mean*janskyConv; cat(paste0("BlanksMeanMean_Jy"," #Mean Value of the Mean Pixel Values in <n> Blank Apertures. In Jansky.\n"))
        newtable[["BlanksMeanSD_Jy"]] = blanks$randMean.SD*janskyConv; cat(paste0("BlanksMeanSD_Jy"," #Standard Deviation of the Mean Pixel Values in <n> Blank Apertures. In Jansky.\n"))
        newtable[["BlanksMeanMAD_Jy"]] = blanks$randMean.MAD*janskyConv; cat(paste0("BlanksMeanMAD_Jy"," #Median Absolute Deviation from Median of the Mean Pixel Values in <n> Blank Apertures. In Jansky.\n"))
        newtable[["BlanksApMean"]] = blanks$randAp.mean; cat(paste0("BlanksApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApSD"]] = blanks$randAp.SD; cat(paste0("BlanksApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApMAD"]] = blanks$randAp.MAD; cat(paste0("BlanksApMAD"," #Median Absolute Deviation from Median of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
      } else {
        newtable[["BlanksMeanMean_units"]] = blanks$randMean.mean; cat(paste0("BlanksMeanMean_units"," #Mean Value of the Mean Pixel Values in <n> Blank Apertures. In image units.\n"))
        newtable[["BlanksMeanSD_units"]] = blanks$randMean.SD; cat(paste0("BlanksMeanSD_units"," #Standard Deviation of the Mean Pixel Values in <n> Blank Apertures. In image units.\n"))
        newtable[["BlanksMeanMAD_units"]] = blanks$randMean.MAD; cat(paste0("BlanksMeanMAD_Jy"," #Median Absolute Deviation from Median of the Mean Pixel Values in <n> Blank Apertures.\n"))
        newtable[["BlanksApMean"]] = blanks$randAp.mean; cat(paste0("BlanksApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApSD"]] = blanks$randAp.SD; cat(paste0("BlanksApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
        newtable[["BlanksApMAD"]] = blanks$randAp.MAD; cat(paste0("BlanksApMAD"," #Median Absolute Deviation from Median of the Integral of the <n> _Individually Masked_ Blank Apertures.\n"))
      }
      newtable[["NumberMeasuredBlanks"]] = blanks$nRand; cat(paste0("NumberMeasuredBlanks"," #Number of Blank Apertures measured per object. Varies object to object due to source masking.\n"))
    }
    if (RanCor) {
      if (Magnitudes) {
        newtable[["RandomsMeanMean_Jy"]] = randoms$randMean.mean*janskyConv; cat(paste0("RandomsMeanMean_Jy"," #Mean Value of the Mean Pixel Values in <n> Random Apertures. In Jansky.\n"))
        newtable[["RandomsMeanSD_Jy"]] = randoms$randMean.SD*janskyConv; cat(paste0("RandomsMeanSD_Jy"," #Standard Deviation of the Mean Pixel Values in <n> Random Apertures. In Jansky.\n"))
        newtable[["RandomsApMean"]] = randoms$randAp.mean; cat(paste0("RandomsApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
        newtable[["RandomsApSD"]] = randoms$randAp.SD; cat(paste0("RandomsApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
      } else {
        newtable[["RandomsMeanMean_units"]] = randoms$randMean.mean; cat(paste0("RandomsMeanMean_units"," #Mean Value of the Mean Pixel Values in <n> Random Apertures. In image units.\n"))
        newtable[["RandomsMeanSD_units"]] = randoms$randMean.SD; cat(paste0("RandomsMeanSD_units"," #Standard Deviation of the Mean Pixel Values in <n> Random Apertures. In image units.\n"))
        newtable[["RandomsApMean"]] = randoms$randAp.mean; cat(paste0("RandomsApMean"," #Mean Value of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
        newtable[["RandomsApSD"]] = randoms$randAp.SD; cat(paste0("RandomsApSD"," #Standard Deviation of the Integral of the <n> _Individually Masked_ Random Apertures.\n"))
      }
    }
    if (!nopsf) {
      newtable[["MinApCorr"]] = ApCorr; cat(paste0("MinApCorr"," #Value of the Multiplicative Minimum Aperture Correction. This correction is important, and can _only_ improve an objects flux determination.\n"))
      newtable[["MaxApCorr"]] = WtCorr; cat(paste0("MaxApCorr"," #Value of the Multiplicative Maximum Aperture Correction. This correction can cause excess flux-boosting systematically as a function of object morphology, and is provided __as a guide only__\n"))
    }
#   Misc.
    newtable[["PhotometryWarning"]] = photWarnFlag; cat(paste0(""," #Some Rudimentary Photometry Warnings. Read Carefully and Take Note!\n"))
    cat(paste0(""," #PhotometryWarnings: Q - Quartered Photometry Warning; 70+% of source flux contained in a single quadrant\n"))
    cat(paste0("","                      S - Sky Estimate Warning; Sky estimate determined in by 3 or less bins\n"))
    cat(paste0("","                      X - Saturation Warning; Pixels in this object aperture are saturated\n"))
    cat(paste0("","                      I - Iteration Warning; This object was lost during iteration because its flux was measured to be <= 0.\n"))
    write.csv(newtable,file=filename, na="-", row.names=FALSE)
  }

  sink(type='output')
  if (!is.null(env)) { detach(env) }
  return=NULL
}

write.deblend.fraction.table <-
function(env=NULL,filename) {
  #Write out a table containing the flux results for each galaxy
  # Load Parameter Space
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  notes<-sub(".csv","_notes.txt",filename)
  sink(type="output",file=notes)
  on.exit(sink(type='output'),add=TRUE)
  today <- Sys.Date()

  cat(paste0("# Notes file for LAMBDAR (v ",packageVersion('LAMBDAR'),") Deblend Fraction Catalogue created ", format(today, format="%B %d %Y"),"\n"))
  cat(paste0("# Below lists all parameters in the accompaning catalogue, and what they represent.\n"))
  if (!psf.weighted) {
    cat(paste0("# NB: This run of LAMBDAR did not request PSF Weighted Apertures. As such, each aperture was passed through a binary-filter which\n",
        "#     converted the contiguous apertures into TopHat apertures. The filter point was requested at ", paste(ap.limit),", which means that\n",
        "#     pixels contained within the ", paste(ap.limit*100),"% contour was set to 1, and the rest were set to 0.\n",
        "# IMPORTANT: Because of this, all references below to the 'Post-Convolution Aperture' (i.e. SFA & DFA) are also after the aperture was binary-filtered.\n"))
  }

##CATALOGUE PARAMETER
  newtable<-data.frame(CATA_INDEX = cat.id)
  colnames(newtable)<-cata.lab; cat(paste0(cata.lab," #Catalogue ID. Duplicates are prepended with 'DuplicatedID_'\n"))
  newtable[[ra.lab]] = cat.ra; cat(paste0(ra.lab," #Right Ascention. Same as provided in the input catalogue\n"))
  newtable[[dec.lab]] = cat.dec; cat(paste0(dec.lab," #Declination. Same as provided in the input catalogue\n"))
  newtable[[theta.lab]] = cat.theta; cat(paste0(theta.lab," #Aperture Orientation Angle. Same as provided in the input catalogue\n"))
  newtable[[semimaj.lab]] = cat.a; cat(paste0(semimaj.lab," #Aperture SemiMajor Axis length, in arcseconds. Same as provided in the input catalogue\n"))
  newtable[[semimin.lab]] = cat.b; cat(paste0(semimin.lab," #Aperture SemiMinor Axis length, in arcseconds. Same as provided in the input catalogue\n"))
  if (exists("contam")) { newtable[[contam.lab]] = contams; cat(paste0(contam.lab," #Contaminant Flag. Same as provided in the input catalogue\n")) }
# IMAGE PARAMETERS
  newtable[["X.imAGE"]] = cat.x; cat(paste0("X.imAGE"," #Aperture Location in Image Coordinates; X-axis\n"))
  newtable[["Y.imAGE"]] = cat.y; cat(paste0("Y.imAGE"," #Aperture Location in Image Coordinates; Y-axis\n"))
  newtable[["StampDiameter_pix"]] = stamplen; cat(paste0("StampDiameter_pix"," #Diameter of the Aperture Stamp, in Pixels.\n"))
# APERTURE SUMS & INTERMEDIATE SUMS
  newtable[["SumPSF"]] = spsf; cat(paste0("SumPSF"," #Integral of the PSF. Should be the same for every object\n"))
  newtable[["SumSA"]] = ssa; cat(paste0("SumSA"," #Integral of the Catalogue Aperture (i.e. the 'Prior') in pixels\n"))
  newtable[["SumSFA"]] = ssfa; cat(paste0("SumSFA"," #Integral of the Post-Convolution Aperture.\n"))
  newtable[["SumDFA"]] = sdfa; cat(paste0("SumDFA"," #Integral of the Deblended Post-Convolution Aperture.\n"))
  newtable[["DeblendFraction"]] = sdfa/ssfa; cat(paste0("DeblendFraction"," #Fractional change in aperture integral caused by deblending the (possibly flux.weighted) apertures.\n"))
  write.csv(newtable,file=filename, na="-", row.names=FALSE)

  sink(type='output')
  if (!is.null(env)) { detach(env) }
  return=NULL
}

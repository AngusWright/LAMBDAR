
#
# Function to get (and plot) the CoG for a single input aperture stamp 
#

get.stamp.cog<-function(outenv=parent.env(environment()), env=NULL) {
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

  #Open Device /*fold*/ {{{
  PlotDev(file=file.path(path.root,path.work,path.out,paste0("COGs/",cat.id[i],".",plot.device)),width=8,height=8,units='in')
  # /*fend*/ }}}
  #Set Layout /*fold*/ {{{
  mar<-par("mar")
  par(mar=mar*c(0.6,0.8,0.5,1))
  layout(rbind(c(1,2),c(3,4),c(0,0)),heights=c(1,1,0.05))
  # /*fend*/ }}}
  #Axes Limits /*fold*/ {{{
  xlims=stamplen[i]*c(-2/3,2/3)*arcsec.per.pix
  ylims=stamplen[i]*c(-2/3,2/3)*arcsec.per.pix
  # /*fend*/ }}}
  #Make Ap and Source Mask Block/Trans matricies /*fold*/ {{{
  apT<-sfa[[i]]
  apT[which(sfa[[i]]==0)]<-NA
  apT[which(sfa[[i]]!=0)]<-1
  apB<-sfa[[i]]
  apB[which(sfa[[i]]==0)]<-1
  apB[which(sfa[[i]]!=0)]<-NA
  if (sourcemask) {
    smB<-sm[mask.stamp.lims[i,1]:mask.stamp.lims[i,2],mask.stamp.lims[i,3]:mask.stamp.lims[i,4]]
    smB[which(smB==0)]<-NA
  } else {
    smB<-1
  }
  # /*fend*/ }}}
  #Plot Image in greyscale /*fold*/ {{{
  Rast<-ifelse(stamplen[i]>100,TRUE,FALSE)
  suppressWarnings(image(x=(seq(1,(diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1))-(x.pix[i]-data.stamp.lims[i,1]))*arcsec.per.pix,y=(seq(1,(diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1))-(y.pix[i]-data.stamp.lims[i,3]))*arcsec.per.pix, z=log10(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=Rast, xlab="", ylab="",xlim=xlims, ylim=ylims,axes=FALSE))
  # /*fend*/ }}}
  #Plot Aperture in Blue /*fold*/ {{{
  suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=(sfa[[i]]), main="Image & Aperture", asp=1, col=hsv(2/3,seq(0,1,length=256)), useRaster=Rast, axes=FALSE, xlab="", ylab="", add=TRUE))
  # /*fend*/ }}}
  #Plot Image in greyscale /*fold*/ {{{
  suppressWarnings(image(x=(seq(1,(diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1))-(x.pix[i]-data.stamp.lims[i,1]))*arcsec.per.pix,y=(seq(1,(diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1))-(y.pix[i]-data.stamp.lims[i,3]))*arcsec.per.pix, z=log10(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=Rast,add=TRUE, xlab="", ylab=""))
  # /*fend*/ }}}
  #Overlay Sourcemask in Green /*fold*/ {{{
  suppressWarnings(image(x=(seq(1,(diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1))-(x.pix[i]-data.stamp.lims[i,1]))*arcsec.per.pix,y=(seq(1,(diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1))-(y.pix[i]-data.stamp.lims[i,3]))*arcsec.per.pix, z=log10(smB*image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]), main="", asp=1, useRaster=Rast,add=TRUE, xlab="", ylab="",col=cm.colors(256)))
  # /*fend*/ }}}
  #Plot +ve flux in aperture in Heat Colours /*fold*/ {{{
  suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=log10(apT*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]), main="", asp=1, col=heat.colors(256), useRaster=Rast,add=TRUE, xlab="", ylab=""))
  # /*fend*/ }}}
  #Plot Sources /*fold*/ {{{
  points(x=(x.pix-x.pix[i]+1)*arcsec.per.pix,y=(y.pix-y.pix[i]+1)*arcsec.per.pix, pch=3)
  # /*fend*/ }}}
  #Label with ID /*fold*/ {{{
  label("topleft",lab=cat.id[i],cex=1.5, col='red')
  # /*fend*/ }}}
  #Label Panel /*fold*/ {{{
  label("topleft",lab="(a)",cex=2.5,inset=c(0.1,0.23))
  # /*fend*/ }}}
  #Draw Axes /*fold*/ {{{
  magaxis(frame.plot=TRUE,main="Image & Aperture",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)",cex.axis=1.2)
  magaxis(side=c(3,4),labels=FALSE)
  # /*fend*/ }}}
  #Generate COGs /*fold*/ {{{
  #Raw Image /*fold*/ {{{
  if (psf.weighted) {
    cog<-get.cog(sfa[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  } else {
    cog<-get.cog(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]],centre=c((diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1)/2, (diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1)/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  }
  # /*fend*/ }}}
  #Sky Subtracted only /*fold*/ {{{
  if (psf.weighted) {
    cog.nosky<-get.cog(sfa[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  } else {
    cog.nosky<-get.cog(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]-skylocal[i],centre=c((diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1)/2, (diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1)/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  }
  # /*fend*/ }}}
  #Deblended only /*fold*/ {{{
  if (psf.weighted) {
    debl.cog<-get.cog(sfa[[i]]*dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  } else {
    debl.cog<-get.cog(dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  }
  # /*fend*/ }}}
  #Deblended and Sky Subtracted /*fold*/ {{{
  if (psf.weighted) {
    debl.cog.nosky<-get.cog(sfa[[i]]*dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  } else {
    debl.cog.nosky<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,poly.degree=Inf,flexible=FALSE)$avg
  }
  # /*fend*/ }}}
  #Deprojected COGs /*fold*/ {{{
  axis.ratio<-cat.b[i]/cat.a[i]
  deprojection.available<-is.finite(axis.ratio) && axis.ratio>0 && axis.ratio<=1
  if (deprojection.available) {
    projection<-c(axis.ratio,theta.offset[i])
    if (psf.weighted) {
      deproj.cog<-get.cog(sfa[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
      deproj.cog.nosky<-get.cog(sfa[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
      deproj.debl.cog<-get.cog(sfa[[i]]*dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
      deproj.debl.cog.nosky<-get.cog(sfa[[i]]*dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
    } else {
      raw.centre<-c((diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1)/2, (diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1)/2)
      deproj.cog<-get.cog(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]],centre=raw.centre,sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
      deproj.cog.nosky<-get.cog(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]-skylocal[i],centre=raw.centre,sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
      deproj.debl.cog<-get.cog(dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
      deproj.debl.cog.nosky<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=projection,poly.degree=Inf,flexible=FALSE)$avg
    }
  }
  # /*fend*/ }}}
  # /*fend*/ }}}
  #Write COGs to text file /*fold*/ {{{
  cog.to.table<-function(x, geometry, deblended, sky.subtracted) {
    data.frame(
      geometry=geometry,
      curve=if (deblended) "deblended" else "blended",
      sky.subtracted=sky.subtracted,
      radius.pix=x$x,
      radius.arcsec=x$x*arcsec.per.pix,
      enclosed.flux=x$y
    )
  }
  cog.table<-rbind(
    cog.to.table(cog, "circular", FALSE, FALSE),
    cog.to.table(debl.cog, "circular", TRUE, FALSE),
    cog.to.table(cog.nosky, "circular", FALSE, TRUE),
    cog.to.table(debl.cog.nosky, "circular", TRUE, TRUE)
  )
  if (deprojection.available) {
    cog.table<-rbind(
      cog.table,
      cog.to.table(deproj.cog, "deprojected", FALSE, FALSE),
      cog.to.table(deproj.debl.cog, "deprojected", TRUE, FALSE),
      cog.to.table(deproj.cog.nosky, "deprojected", FALSE, TRUE),
      cog.to.table(deproj.debl.cog.nosky, "deprojected", TRUE, TRUE)
    )
  }
  write.table(
    cog.table,
    file=file.path(path.root,path.work,path.out,paste0("COGs/",cat.id[i],".txt")),
    sep="\t",
    row.names=FALSE,
    quote=FALSE,
    na="NA"
  )
  # /*fend*/ }}}
  #Plot COGs /*fold*/ {{{
  if (magnitudes) {
    #Plot in Magnitude Space /*fold*/ {{{
    #Get y limits /*fold*/ {{{
    suppressWarnings(ylim<-c((-2.5*(log10(dfaflux[i])-log10(ab.vega.flux))+mag.zp)+3, (-2.5*(log10(dfaflux[i])-log10(ab.vega.flux))+mag.zp)-3))
    # /*fend*/ }}}
    #If a limit in ±Inf, plot around median /*fold*/ {{{
    if (!all(is.finite(ylim))) { ylim=median(-2.5*(log10(cog$y)-log10(ab.vega.flux))+mag.zp, na.rm=TRUE)+c(-1,3) }
    if (!all(is.finite(ylim))) { warning("Cog flux is always < 0; Using arbitrary plot limits"); ylim=18+c(-1,3) }
    # /*fend*/ }}}
    #Plot Raw COG /*fold*/ {{{
    suppressWarnings(magplot(x=cog$x*arcsec.per.pix, y=-2.5*(log10(cog$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Magnitude",
            ylim=ylim,xlim=c(0,max(xlims)),type='l',lty=2,main="Curve of Growth",cex.axis=1.2))
    magaxis(side=c(3,4),labels=FALSE)
    # /*fend*/ }}}
    #Add Lines showing fluxes /*fold*/ {{{
    #Undeblended /*fold*/ {{{
    suppressWarnings(abline(h=-2.5*(log10(sfaflux[i])-log10(ab.vega.flux))+mag.zp, lwd=1, col='orange', lty=1))
    # /*fend*/ }}}
    #Deblended /*fold*/ {{{
    suppressWarnings(abline(h=-2.5*(log10(dfaflux[i])-log10(ab.vega.flux))+mag.zp, lwd=1, col='green', lty=1))
    # /*fend*/ }}}
    #Redraw Raw Cog /*fold*/ {{{
    suppressWarnings(lines(x=cog$x*arcsec.per.pix, y=-2.5*(log10(cog$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='grey',lty=2))
    # /*fend*/ }}}
    #Draw Deblended Cog /*fold*/ {{{
    suppressWarnings(lines(x=debl.cog$x*arcsec.per.pix, y=-2.5*(log10(debl.cog$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='black',lty=2))
    # /*fend*/ }}}
    #Draw Sky Subtracted Cog /*fold*/ {{{
    suppressWarnings(lines(x=cog.nosky$x*arcsec.per.pix, y=-2.5*(log10(cog.nosky$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='grey',lty=1))
    # /*fend*/ }}}
    #Draw Sky Subtracted & Deblended Cog /*fold*/ {{{
    suppressWarnings(lines(x=debl.cog.nosky$x*arcsec.per.pix, y=-2.5*(log10(debl.cog.nosky$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='black',lty=1))
    # /*fend*/ }}}
    #Draw Legend /*fold*/ {{{
    legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended ApMag","Deblended ApMag"),lty=c(2,2,1,1,1,1),
           col=c('grey','black','grey','black','orange','green'),pch=-1, cex=1.2)
    # /*fend*/ }}}
    #Note Half Light Radius /*fold*/ {{{
    if (do.sky.est) {
      hlr<-which(debl.cog.nosky$y>=max(debl.cog.nosky$y,na.rm=TRUE)/2)
      if (deprojection.available) {
        dhlr<-which(deproj.debl.cog.nosky$y>max(debl.cog.nosky$y,na.rm=TRUE)/2)
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog.nosky$x[hlr]),digits=2),"\nDeprojected:",round(min(deproj.debl.cog.nosky$x[dhlr]),digits=2)),cex=1.2)
      } else {
        dhlr<-hlr
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog.nosky$x[hlr]),digits=2)),cex=1.2)
      }
    } else {
      hlr<-which(debl.cog$y>=max(debl.cog$y,na.rm=TRUE)/2)
      if (deprojection.available) {
        dhlr<-which(deproj.debl.cog$y>max(debl.cog$y,na.rm=TRUE)/2)
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog$x[hlr]),digits=2),"\nDeprojected:",round(min(deproj.debl.cog$x[dhlr]),digits=2)),cex=1.2)
      } else {
        dhlr<-hlr
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog$x[hlr]),digits=2)),cex=1.2)
      }
    }
    # /*fend*/ }}}
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else {
    #Plot in Flux Space /*fold*/ {{{
    #Plot Raw Cog /*fold*/ {{{
    ylim<-range(cog$y)
    #If a limit in ±Inf, plot around median /*fold*/ {{{
    if (!all(is.finite(ylim))) { ylim=median(cog$y, na.rm=TRUE)+c(-1E2,1E2) }
    magplot(x=cog$x*arcsec.per.pix, y=cog$y, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Flux",ylim=ylim,xlim=c(0,max(xlims)),main="Curve of Growth",type='l',cex.axis=1.2)
    magaxis(side=c(3,4),labels=FALSE)
    # /*fend*/ }}}
    # /*fend*/ }}}
    #Draw Lines showing fluxes /*fold*/ {{{
    abline(h=dfaflux[i], lwd=1, col='green')
    abline(h=sfaflux[i], lwd=1, col='orange', lty=2)
    # /*fend*/ }}}
    #Redraw Raw Cog /*fold*/ {{{
    lines(x=cog$x*arcsec.per.pix, y=cog$y, pch=20, col='grey',lty=2)
    # /*fend*/ }}}
    #Draw Deblended Cog /*fold*/ {{{
    lines(x=debl.cog$x*arcsec.per.pix, y=debl.cog$y, pch=20, col='black',lty=2)
    # /*fend*/ }}}
    #Draw Sky Subtracted Cog /*fold*/ {{{
    lines(x=cog.nosky$x*arcsec.per.pix, y=cog.nosky$y, pch=20, col='grey',lty=1)
    # /*fend*/ }}}
    #Draw Sky Subtracted & Deblended Cog /*fold*/ {{{
    lines(x=debl.cog.nosky$x*arcsec.per.pix, y=debl.cog.nosky$y, pch=20, col='black',lty=1)
    # /*fend*/ }}}
    legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended Flux","Deblended Flux"),lty=c(2,2,1,1,1,1),
           col=c('grey','black','grey','black','orange','green'),pch=-1, cex=1.2)
    #Note Half Light Radius /*fold*/ {{{
    if (do.sky.est) {
      hlr<-which(debl.cog.nosky$y>=max(cog$y,na.rm=TRUE)/2)
      if (deprojection.available) {
        dhlr<-which(deproj.debl.cog.nosky$y>max(debl.cog$y,na.rm=TRUE)/2)
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog.nosky$x[hlr]),digits=2),"\nDeprojected:",round(min(deproj.debl.cog.nosky$x[dhlr]),digits=2)))
      } else {
        dhlr<-hlr
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog.nosky$x[hlr]),digits=2)))
      }
    } else {
      hlr<-which(debl.cog$y>=max(cog$y,na.rm=TRUE)/2)
      if (deprojection.available) {
        dhlr<-which(deproj.debl.cog$y>max(debl.cog$y,na.rm=TRUE)/2)
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog$x[hlr]),digits=2),"\nDeprojected:",round(min(deproj.debl.cog$x[dhlr]),digits=2)))
      } else {
        dhlr<-hlr
        label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog$x[hlr]),digits=2)))
      }
    }
    # /*fend*/ }}}
    # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #Label Panel /*fold*/ {{{
  label("topleft",lab="(b)",cex=2.5,inset=c(0.1,0.23))
  # /*fend*/ }}}
  #Plot the Deblended Image /*fold*/ {{{
  nc<-length(ap.lims.data.map[i,1]:ap.lims.data.map[i,2])
  nr<-length(ap.lims.data.map[i,3]:ap.lims.data.map[i,4])
  Rast<-ifelse(stamplen[i]>100,TRUE,FALSE)
  z<-dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i])
  if (any(is.finite(z))) { 
    suppressWarnings(z<-matrix(magmap(z,stretch='asinh')$map,ncol=nc,nrow=nr))
  }
  image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apB, main="Image x Weight Matrix", asp=1, col=grey.colors(256), useRaster=Rast, xlab="", ylab="", axes=FALSE, xlim=xlims, ylim=ylims)
  # /*fend*/ }}}
  #Overlay the Aperture /*fold*/ {{{
  suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apT, main="Image x Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=Rast, xlab="", ylab="", axes=FALSE, xlim=xlims, ylim=ylims,add=TRUE))
  # /*fend*/ }}}
  #Draw the projected half-light ellipse {{{
  if (deprojection.available) {
    deproj.hlr.cog<-if (do.sky.est) deproj.debl.cog.nosky else deproj.debl.cog
    lines(ellipse(a=min(deproj.hlr.cog$x[dhlr])*arcsec.per.pix,e=1-axis.ratio,pa=90-theta.offset[i],x0=(cat.x[i]%%1)*arcsec.per.pix,y0=(cat.y[i]%%1)*arcsec.per.pix),col='black',lty=3,lwd=2)
  }
  #}}}
  #Draw the Axes and scalebar /*fold*/ {{{
  magaxis(frame.plot=TRUE,main="Image x Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)",cex.axis=1.2)
  magaxis(side=c(3,4),labels=FALSE)
  # /*fend*/ }}}
  #Label Panel /*fold*/ {{{
  label("topleft",lab="(c)",cex=2.5,inset=c(0.1,0.23))
  # /*fend*/ }}}
  #Plot the Deblend Matrix /*fold*/ {{{
  z=dbw[[i]]
  image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apB, main="Weight Matrix", asp=1, col=grey.colors(256), useRaster=Rast, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=ylims)
  # /*fend*/ }}}
  #Overlay the Aperture /*fold*/ {{{
  suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apT, main="Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=Rast, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=ylims,add=TRUE))
  # /*fend*/ }}}
  #Plot Sources /*fold*/ {{{
  points(x=(x.pix-x.pix[i]+1)*arcsec.per.pix,y=(y.pix-y.pix[i]+1)*arcsec.per.pix, pch=3,col=hsv(v=0,a=0.3))
  # /*fend*/ }}}
  #Draw the Axes and scalebar /*fold*/ {{{
  magaxis(frame.plot=TRUE,main="Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)",cex.axis=1.2)
  magaxis(side=c(3,4),labels=FALSE)
  # /*fend*/ }}}
  #Label Panel /*fold*/ {{{
  label("topleft",lab="(d)",cex=2.5,inset=c(0.1,0.23))
  # /*fend*/ }}}
  #Close the file /*fold*/ {{{
  if (!grepl('x11',plot.device,ignore.case=TRUE)) { dev.off() }
  # /*fend*/ }}}
  #Plot deprojected COGs /*fold*/ {{{
  if (deprojection.available) {
    PlotDev(file=file.path(path.root,path.work,path.out,paste0("COGs/",cat.id[i],"_deprojected.",plot.device)),width=8,height=8,units='in')
    par(mar=mar)
    deproj.xlim<-c(0,max(c(deproj.cog$x,deproj.debl.cog$x,deproj.cog.nosky$x,deproj.debl.cog.nosky$x),na.rm=TRUE)*arcsec.per.pix)
    if (magnitudes) {
      to.mag<-function(x) -2.5*(log10(x)-log10(ab.vega.flux))+mag.zp
      suppressWarnings(deproj.ylim<-c(to.mag(dfaflux[i])+3,to.mag(dfaflux[i])-3))
      if (!all(is.finite(deproj.ylim))) {
        suppressWarnings(deproj.ylim<-median(to.mag(deproj.cog$y),na.rm=TRUE)+c(-1,3))
      }
      if (!all(is.finite(deproj.ylim))) { deproj.ylim<-18+c(-1,3) }
      suppressWarnings(magplot(deproj.cog$x*arcsec.per.pix,to.mag(deproj.cog$y),type='l',lty=2,col='grey',xlab="Deprojected semi-major radius (arcsec)",ylab="Enclosed Magnitude",xlim=deproj.xlim,ylim=deproj.ylim,main="Deprojected Curve of Growth",cex.axis=1.2))
      suppressWarnings(lines(deproj.debl.cog$x*arcsec.per.pix,to.mag(deproj.debl.cog$y),lty=2,col='black'))
      suppressWarnings(lines(deproj.cog.nosky$x*arcsec.per.pix,to.mag(deproj.cog.nosky$y),lty=1,col='grey'))
      suppressWarnings(lines(deproj.debl.cog.nosky$x*arcsec.per.pix,to.mag(deproj.debl.cog.nosky$y),lty=1,col='black'))
      suppressWarnings(abline(h=to.mag(sfaflux[i]),lwd=1,col='orange'))
      suppressWarnings(abline(h=to.mag(dfaflux[i]),lwd=1,col='green'))
      legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended ApMag","Deblended ApMag"),lty=c(2,2,1,1,1,1),col=c('grey','black','grey','black','orange','green'),pch=-1,cex=1.2)
    } else {
      deproj.ylim<-range(c(deproj.cog$y,deproj.debl.cog$y,deproj.cog.nosky$y,deproj.debl.cog.nosky$y),finite=TRUE)
      magplot(deproj.cog$x*arcsec.per.pix,deproj.cog$y,type='l',lty=2,col='grey',xlab="Deprojected semi-major radius (arcsec)",ylab="Enclosed Flux",xlim=deproj.xlim,ylim=deproj.ylim,main="Deprojected Curve of Growth",cex.axis=1.2)
      lines(deproj.debl.cog$x*arcsec.per.pix,deproj.debl.cog$y,lty=2,col='black')
      lines(deproj.cog.nosky$x*arcsec.per.pix,deproj.cog.nosky$y,lty=1,col='grey')
      lines(deproj.debl.cog.nosky$x*arcsec.per.pix,deproj.debl.cog.nosky$y,lty=1,col='black')
      abline(h=sfaflux[i],lwd=1,col='orange')
      abline(h=dfaflux[i],lwd=1,col='green')
      legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended Flux","Deblended Flux"),lty=c(2,2,1,1,1,1),col=c('grey','black','grey','black','orange','green'),pch=-1,cex=1.2)
    }
    mtext(paste0("Axis ratio b/a = ",round(axis.ratio,3),"; PA = ",round(theta.offset[i],2)," deg"),side=3,line=0.25,cex=0.9)
    magaxis(side=c(3,4),labels=FALSE)
    if (!grepl('x11',plot.device,ignore.case=TRUE)) { dev.off() }
  }
  # /*fend*/ }}}
}

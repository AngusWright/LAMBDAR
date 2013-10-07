make_rgb <-
function (struct=NA, file=NA, points=FALSE, cat=NA, lo=0.4, hi=0.995, log=TRUE, bad=0, flip=FALSE, outpng=NA, r=3, g=2, b=1, title="DEFAULT", ...) {

   if (all(is.na(struct))&(is.na(file))) { stop("No File or Image List Provided") }
   if (!is.na(file)) { struct<-read.fits(file, hdu=0) }
   if (!is.na(outpng)) { png(file=outpng, width=max(1200, length(struct$dat[[r]][1,])+100), height=max(1200, length(struct$dat[[r]][1,])+100)) }
   plot.new()
   if (!(all(is.nan(struct$dat[[r]]))|all(is.nan(struct$dat[[g]]))|all(is.nan(struct$dat[[b]])))){
     cols<-rgb(magmap(struct$dat[[r]],bad=bad,lo=lo,hi=hi,log=log,flip=flip,...)$map, 
               magmap(struct$dat[[g]],bad=bad,lo=lo,hi=hi,log=log,flip=flip,...)$map, 
               magmap(struct$dat[[b]],bad=bad,lo=lo,hi=hi,log=log,flip=flip,...)$map) 
     plot.window(xlim=range(struct$x[[r]],na.rm=TRUE), 
                 ylim=range(struct$y[[r]],na.rm=TRUE), 
                  asp=1/cos(median(struct$y[[r]], na.rm=TRUE)*pi/180))
     image(x=seq(min(struct$x[[r]]),max(struct$x[[r]]),len=length(struct$x[[r]])),
           y=seq(min(struct$y[[r]]),max(struct$y[[r]]),len=length(struct$y[[r]])),
           z=matrix(1:length(struct$dat[[r]]), ncol=dim(struct$dat[[r]])[2], byrow=T), 
         col=cols, axes=FALSE, add=TRUE, useRaster=TRUE)
     #print(par("usr"))
     plot.window(xlim=c(struct$x[[r]][1],struct$x[[r]][length(struct$x[[r]])]), 
                 ylim=c(struct$y[[r]][1],struct$y[[r]][length(struct$y[[r]])]), 
                  asp=1/cos(median(struct$y[[r]], na.rm=TRUE)*pi/180))
     magaxis(xlab="RA/deg",ylab="Dec/deg")
     if (points) {
       if (is.na(cat)) { 
         warning("No CSV Catalogue Supplied with points flag") 
       } else {
         try.load=try(name<-load(cat))
         if (class(try.load)=='try error') {
           warning("Failed to read points catalogue")
         } else {
           cat<-get(name)
           if (is.null(cat$PLOT_MAIN)) {
             if (!flip) { points(cat$RA, cat$DEC, cex=magmap(cat$R_PETRO,lo=15,hi=20,range=c(7,3),type='num')$map, col='lightgreen') }
             if ( flip) { points(cat$RA, cat$DEC, cex=magmap(cat$R_PETRO,lo=15,hi=20,range=c(7,3),type='num')$map, col='red'       ) }
           } else if (is.null(cat$PLOT_REDSHIFT)){
             if (!flip) { points(cat$RA[which( cat$PLOT_MAIN)], cat$DEC[which( cat$PLOT_MAIN)], cex=magmap(cat$R_PETRO[which( cat$PLOT_MAIN)],lo=15,hi=20,range=c(7,3),type='num')$map,        col='lightgreen', lwd=3)
                          points(cat$RA[which(!cat$PLOT_MAIN)], cat$DEC[which(!cat$PLOT_MAIN)], cex=magmap(cat$R_PETRO[which(!cat$PLOT_MAIN)],lo=15,hi=20,range=c(7,3),type='num')$map, pch=4, col='lightgreen', lwd=3) }
             if ( flip) { points(cat$RA[which( cat$PLOT_MAIN)], cat$DEC[which( cat$PLOT_MAIN)], cex=magmap(cat$R_PETRO[which( cat$PLOT_MAIN)],lo=15,hi=20,range=c(7,3),type='num')$map,        col='red'       , lwd=3) 
                          points(cat$RA[which(!cat$PLOT_MAIN)], cat$DEC[which(!cat$PLOT_MAIN)], cex=magmap(cat$R_PETRO[which(!cat$PLOT_MAIN)],lo=15,hi=20,range=c(7,3),type='num')$map, pch=4, col='red'       , lwd=3) }
           } else {
             if (!flip) { points(cat$RA[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))], cat$DEC[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))], 
                 cex=magmap(cat$R_PETRO[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))],lo=15,hi=20,range=c(7,3),type='num')$map,        col='lightgreen', lwd=3)
                          points(cat$RA[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))], cat$DEC[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))],
                 cex=magmap(cat$R_PETRO[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))],lo=15,hi=20,range=c(7,3),type='num')$map,        col='orange', lwd=3)
                          points(cat$RA[which(!(cat$PLOT_MAIN))], cat$DEC[which(!(cat$PLOT_MAIN))],
                 cex=magmap(cat$R_PETRO[which(!(cat$PLOT_MAIN))],lo=15,hi=20,range=c(7,3),type='num')$map, pch=4, col='red', lwd=3) }
             if ( flip) { points(cat$RA[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))], cat$DEC[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))],
                 cex=magmap(cat$R_PETRO[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))],lo=15,hi=20,range=c(7,3),type='num')$map,        col='darkgreen'       , lwd=3) 
                          points(cat$RA[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))], cat$DEC[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))],
                 cex=magmap(cat$R_PETRO[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))],lo=15,hi=20,range=c(7,3),type='num')$map,        col='darkorange'       , lwd=3) 
                          points(cat$RA[which(!(cat$PLOT_MAIN))], cat$DEC[which(!(cat$PLOT_MAIN))],
                 cex=magmap(cat$R_PETRO[which(!(cat$PLOT_MAIN))],lo=15,hi=20,range=c(7,3),type='num')$map, pch=4, col='red'       , lwd=3) }
           }
         }
       }
     } 
     plot.window(xlim=c(0,1), ylim=c(0,1), asp=1)
     if (flip) { text(0.15, 0.9, title, cex=3, col="darkred") }
     else { text(0.15, 0.9, title, cex=3, col="green") }
     if (!is.na(outpng)) { dev.off() ; return=NULL }
   } else {
     plot.new()
     image(array(1,dim=dim(struct$dat[[r]])), axes=FALSE,  col=grey(0.5))
     plot.window(xlim=c(0,1), ylim=c(0,1), asp=1)
     if (flip) { text(0.15, 0.9, title, cex=3, col="darkred") 
       text(0.4, 0.5, "MISSING DATA", cex=5, col="darkred")
     } else { text(0.15, 0.9, title, cex=3, col="green")
       text(0.4, 0.5, "MISSING DATA", cex=5, col="red")
     }
     plot.window(xlim=c(struct$x[[r]][1],struct$x[[r]][length(struct$x[[r]])]), 
                 ylim=c(struct$y[[r]][1],struct$y[[r]][length(struct$y[[r]])]), 
                  asp=1/cos(median(struct$y[[r]], na.rm=TRUE)*pi/180))
     magaxis(xlab="RA/deg",ylab="Dec/deg")
     if (!any(is.na(outpng))) { dev.off() } 
   }  
   return=list(vals=matrix(1:length(struct$dat[[1]]), ncol=dim(struct$dat[[1]])[2], byrow=T), cols=cols)
}

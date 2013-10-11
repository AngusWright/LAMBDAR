make_greyscale <-
function (struct=NA, file=NA, points=FALSE, cat=NA, lo=0.1, hi=0.998, log=TRUE, bad=0, flip=FALSE, outpng=NA, title="",...) {

   if (all(is.na(struct))&(is.na(file))) { stop("No File or Image List Provided") }
   if (!is.na(file)) { struct<-read.fits(file, hdu=0)}
   if (points) {
     if (is.na(cat)) {
       warning("No CSV Catalogue Supplied with points flag")
       points=FALSE
     } else {
       try.read=try(catName<-load(cat))
       if (class(try.read)=='try error') {
         warning("Failed to read points catalogue")
         points=FALSE
       }
     }
   }
   for (i in 1:length(struct$dat) ) {
     if (!any(is.na(outpng))) { png(filename=outpng[i]) }
     plot.new()
     if (!all(is.nan(struct$dat[[i]]))) {
       map<-magmap(struct$dat[[i]],bad=bad,lo=lo,hi=hi,log=log,flip=flip,...)$map
       #Rotate matrix to correct x/y axis orientation
       map<-matrix(map, ncol=dim(map)[2], byrow=T)
       col<-rgb(map,map,map)
       plot.window(xlim=range(struct$x[[i]],na.rm=TRUE),
                   ylim=range(struct$y[[i]],na.rm=TRUE),
                    asp=1/cos(median(struct$y[[i]], na.rm=TRUE)*pi/180))
       image(x=seq(min(struct$x[[i]]),max(struct$x[[i]]),len=length(struct$x[[i]])),
             y=seq(min(struct$y[[i]]),max(struct$y[[i]]),len=length(struct$y[[i]])),
             #z=matrix(1:length(struct$dat[[i]]), ncol=dim(struct$dat[[i]])[2], byrow=T),col=col, axes=FALSE, add=TRUE, useRaster=TRUE)
             z=map,col=grey((0:128)/128), axes=FALSE, add=TRUE, useRaster=TRUE, main=title)
       #print(par("usr"))
       plot.window(xlim=c(struct$x[[i]][1],struct$x[[i]][length(struct$x[[i]])]),
                   ylim=c(struct$y[[i]][1],struct$y[[i]][length(struct$y[[i]])]),
                    asp=1/cos(median(struct$y[[i]], na.rm=TRUE)*pi/180))
       magaxis(xlab="RA/deg",ylab="Dec/deg", cex=2)
       if (points) {
         cat<-get(catName)
           if (is.null(cat$PLOT_MAIN)) {
             if (!flip) { points(cat$RA, cat$DEC, cex=magmap(cat$R_PETRO,lo=15,hi=20,range=c(7,3),type='num')$map, col='lightgreen',lwd=3) }
             if ( flip) { points(cat$RA, cat$DEC, cex=magmap(cat$R_PETRO,lo=15,hi=20,range=c(7,3),type='num')$map, col='red'       ,lwd=3) }
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
                 cex=magmap(cat$R_PETRO[which(!(cat$PLOT_MAIN))],lo=15,hi=20,range=c(7,3),type='num')$map, pch=4, col='pink', lwd=3) }
             if ( flip) { points(cat$RA[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))], cat$DEC[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))],
                 cex=magmap(cat$R_PETRO[which( (cat$PLOT_MAIN)&(cat$PLOT_REDSHIFT))],lo=15,hi=20,range=c(7,3),type='num')$map,        col='darkgreen'       , lwd=3)
                          points(cat$RA[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))], cat$DEC[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))],
                 cex=magmap(cat$R_PETRO[which( (cat$PLOT_MAIN)&(!cat$PLOT_REDSHIFT))],lo=15,hi=20,range=c(7,3),type='num')$map,        col='darkorange'       , lwd=3)
                          points(cat$RA[which(!(cat$PLOT_MAIN))], cat$DEC[which(!(cat$PLOT_MAIN))],
                 cex=magmap(cat$R_PETRO[which(!(cat$PLOT_MAIN))],lo=15,hi=20,range=c(7,3),type='num')$map, pch=4, col='red'       , lwd=3) }
           }
       }
       plot.window(xlim=c(0,1), ylim=c(0,1), asp=1)
       if (flip) { text(0.15, 0.9, title, cex=3, col="darkred") }
       else { text(0.15, 0.9, title, cex=3, col="green") }
       if (!any(is.na(outpng))) { dev.off() }
     } else {
       plot.new()
       image(array(1,dim=dim(struct$dat[[i]])), axes=FALSE,  col=grey(0.5))
       plot.window(xlim=c(0,1), ylim=c(0,1), asp=1)
       if (flip) { text(0.15, 0.9, title, cex=3, col="darkred")
       text(0.4, 0.5, "NO DATA", cex=5, col="darkred")
       } else { text(0.15, 0.9, title, cex=3, col="green")
       text(0.4, 0.5, "NO DATA", cex=5, col="red")
       }
       plot.window(xlim=c(struct$x[[i]][1],struct$x[[i]][length(struct$x[[i]])]),
                   ylim=c(struct$y[[i]][1],struct$y[[i]][length(struct$y[[i]])]),
                    asp=1/cos(median(struct$y[[i]], na.rm=TRUE)*pi/180))
       magaxis(xlab="RA/deg",ylab="Dec/deg")
       if (!any(is.na(outpng))) { dev.off() }
     }
   }
   return=NULL
}

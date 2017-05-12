fill.grid<-
function(ra, dec, buffer=0.1, rad=sqrt(3)/2, buffer.edge=TRUE, plot=TRUE, main="Block Coverage Map") {
  if (buffer.edge) {
    ra[which.min(ra)]<-min(ra)-buffer
    ra[which.max(ra)]<-max(ra)+buffer
    dec[which.min(dec)]<-min(dec)-buffer
    dec[which.max(dec)]<-max(dec)+buffer
  }
  ra<-range(ra)
  dec<-range(dec)
  if (rad<=0) {
    warning("Bad Radius supplied. Using sqrt(3)/2")
    rad=sqrt(3)/2
  }
  if (buffer > rad) {
    warning("Buffer cannot be greater than Radius, as the algorithm will be non-convergent...")
    buffer<-rad*0.75
  }
  #Number of RA blocks {{{
  nblocks<-ceiling(diff(ra)/(rad*2))
  limit<-3*nblocks
  finished=FALSE
  while (!finished) {
    if (diff(ra)>nblocks*(rad*2)-(nblocks-1)*(buffer)) {
      nblocks<-nblocks+1
    } else {
      finished=TRUE
    }
  }
  nblocks.ra<-nblocks
  #}}}
  #Number of Dec Blocks {{{
  nblocks<-ceiling(diff(dec)/(rad*2))
  finished=FALSE
  while (!finished) {
    if (diff(dec)>nblocks*(rad*2)-(nblocks-1)*(buffer)) {
      nblocks<-nblocks+1
    } else {
      finished=TRUE
    }
  }
  nblocks.dec<-nblocks
  #}}}
  ra_step=NULL
  for (i in 1:nblocks.ra) {
    ra_step<-c(ra_step,min(ra)+(i-0.5)*(rad*2)-(i-1)*buffer)
  }
  dec_step=NULL
  for (i in 1:nblocks.dec) {
    dec_step<-c(dec_step,min(dec)+(i-0.5)*(rad*2)-(i-1)*buffer)
  }
  overflow<-max(ra_step)+rad-max(ra)
  ra_step=ra_step-overflow/2
  overflow<-max(dec_step)+rad-max(dec)
  dec_step=dec_step-overflow/2
  grid=expand.grid(ra_step, dec_step)
  if (plot) {
   plot(expand.grid(ra,dec)[c(1,2,4,3,1),c(1,2)], type='l', asp=1, xlab="RA", ylab="DEC", main=main, xlim=range(ra_step)+c(-rad,rad), ylim=range(dec_step)+c(-rad,rad))
   for (i in 1:length(grid[,1])) {
     rect(xleft=grid[i,1]-rad, ybottom=grid[i,2]-rad, xright=grid[i,1]+rad, ytop=grid[i,2]+rad, lty=3)
   }
   points(grid, pch=3)
  label('bottomleft',(paste("Number of RA  Blocks:", nblocks.ra, "\nNumber of DEC Blocks:", nblocks.dec)))
  }
  return=grid
}

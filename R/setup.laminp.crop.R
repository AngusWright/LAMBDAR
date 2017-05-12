#
#
#
#
#
#

setup.laminp.crop<-function(laminp.imagelist,laminp.workdirlist,crop.rad,buffer,other.laminp.filelists) { 

  #Initialise the LAMINP_ files
  system("rm -f LAMINP_CROPRA.dat LAMINP_CROPDEC.dat LAMINP_CROPRAD.dat")
  system("touch LAMINP_CROPRA.dat LAMINP_CROPDEC.dat LAMINP_CROPRAD.dat")
  #If other LAMINP filelist, move them to a .bak extension 
  if (!missing(other.laminp.filelists)) {
    system(paste0('mv ',other.laminp.filelists,' ',other.laminp.filelists,".bak",collapse=' ; '))
    system(paste0('touch ',other.laminp.filelists,collapse=' ; '))
  }
  images<-read.table(file=laminp.imagelist,header=FALSE,stringsAsFactors=FALSE)$V1
  system(paste0('mv ',laminp.imagelist,' ',laminp.imagelist,".bak",collapse=' ; '))
  system(paste0('touch ',laminp.imagelist,collapse=' ; '))
  #Working directory list
  if (missing(laminp.workdirlist)) {
    workdirlist<-rep("./",length(images))
  } else {
    workdirlist<-read.table(file=laminp.workdirlist,header=FALSE,stringsAsFactors=FALSE)$V1
    system(paste0('mv ',laminp.workdirlist,' ',laminp.workdirlist,".bak",collapse=' ; '))
    system(paste0('touch ',laminp.workdirlist,collapse=' ; '))
    if (length(workdirlist)!=length(images)) { 
      stop("workdir list and image list are of different lengths!")
    }
  }
  #Loop over all images
  for (i in 1:length(images)) {
    #get the image RA and DEC limits 
    hdr<-read.astrometry(paste0(workdirlist[i],'/',images[i]))
    lims<-xy.to.ad(c(1,1,hdr$NAXIS[1],hdr$NAXIS[1]),c(1,hdr$NAXIS[2],hdr$NAXIS[2],1),hdr)
    #calculate the grid of crops
    grid<-fill.grid(ra=lims[,1],dec=lims[,2],rad=crop.rad,buffer=buffer,plot=FALSE)
    #Add to LAMINP_ lists
    write.table(file='LAMINP_CROPRA.dat' ,grid[,1],quote=FALSE,row.names=FALSE,append=TRUE,col.names=FALSE)
    write.table(file='LAMINP_CROPDEC.dat',grid[,2],quote=FALSE,row.names=FALSE,append=TRUE,col.names=FALSE)
    write.table(file='LAMINP_CROPRAD.dat',rep(crop.rad,length(grid[,1])),quote=FALSE,row.names=FALSE,append=TRUE,col.names=FALSE)
    write.table(file=laminp.imagelist,rep(images[i],length(grid[,1])),quote=FALSE,row.names=FALSE,append=TRUE,col.names=FALSE)
    #Update LAMINP_workdirs file
    if (!missing(laminp.workdirlist)) { 
      write.table(file=laminp.workdirlist,rep(workdirlist[i],length(grid[,1])),quote=FALSE,row.names=FALSE,append=TRUE,col.names=FALSE)
    }
    #Update other LAMINP_ files 
    for (j in 1:length(other.laminp.filelists)) { 
      val<-read.table(file=paste0(other.laminp.filelists[j],'.bak'),header=FALSE,stringsAsFactors=FALSE)$V1[i]
      write.table(file=other.laminp.filelists[j],rep(val,length(grid[,1])),quote=FALSE,row.names=FALSE,append=TRUE,col.names=FALSE)
    } 
  }
}

iterapint <-
function(x,y,xstep,ystep,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=T,upres=2,itersteps=9,pixscale=FALSE,peakscale=FALSE,peak=1){

if(majax>0){

 if(itersteps>0){

  origx=x
  origy=y
  origxstep=xstep
  origystep=ystep

  mastercat={}
  weight=1
  newID=1:length(x)

  for(i in 1:itersteps){
   check=checkgrid(x=x,y=y,xstep=xstep,ystep=ystep,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
   if(any(check$logic[,'empty'])){mastercat=rbind(mastercat,cbind(newID[check$logic[,'empty']],0))}
   if(any(check$logic[,'full'])){mastercat=rbind(mastercat,cbind(newID[check$logic[,'full']],weight))} 
   regrid=reexpand.grid(x[check$logic[,'part']],y[check$logic[,'part']],xstep=xstep,ystep=ystep,id=newID[check$logic[,'part']],upres=upres)
   x=regrid$x
   y=regrid$y
   newID=regrid$id
   xstep=regrid$xstep
   ystep=regrid$ystep
   weight=weight/(upres^2)
  }

  check=checkgrid(x=x,y=y,xstep=xstep,ystep=ystep,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
  if(any(check$logic[,'empty'])){mastercat=rbind(mastercat,cbind(newID[check$logic[,'empty']],0))}
  if(any(check$logic[,'full'])){mastercat=rbind(mastercat,cbind(newID[check$logic[,'full']],weight))}
  if(any(check$logic[,'mid']==F & check$logic[,'part'])){mastercat=rbind(mastercat,cbind(newID[check$logic[,'mid']==F & check$logic[,'part']],0))}
  if(any(check$logic[,'mid'] & check$logic[,'part'])){mastercat=rbind(mastercat,cbind(newID[check$logic[,'mid'] & check$logic[,'part']],weight))}

  rebinID=as.numeric(names(table(mastercat[,1])[table(mastercat[,1])>1]))

  sumcat={}

  for(i in 1:length(rebinID)){
   sumcat=rbind(sumcat,cbind(rebinID[i],sum(mastercat[mastercat[,1]==rebinID[i],2])))
  }
  mastercat=rbind(mastercat[!mastercat[,1] %in% rebinID,],sumcat)
  mastercat=mastercat[order(mastercat[,1]),2]
  mastercat=cbind(origx,origy,mastercat)
  if(pixscale){mastercat[,3]=mastercat[,3]*origxstep*origystep}
 }else{
  mastercat={}
  weight=1

  check=checkgrid(x=x,y=y,xstep=xstep,ystep=ystep,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
  if(any(check$logic[,'empty'])){mastercat=rbind(mastercat,cbind(which(check$logic[,'empty']),0))}
  if(any(check$logic[,'full'])){mastercat=rbind(mastercat,cbind(which(check$logic[,'full']),weight))}
  if(any(check$logic[,'mid']==F & check$logic[,'part'])){mastercat=rbind(mastercat,cbind(which(check$logic[,'mid']==F & check$logic[,'part']),0))}
  if(any(check$logic[,'mid'] & check$logic[,'part'])){mastercat=rbind(mastercat,cbind(which(check$logic[,'mid'] & check$logic[,'part']),weight))}

  mastercat=mastercat[order(mastercat[,1]),2]
  mastercat=cbind(x,y,mastercat)
  if(pixscale){mastercat[,3]=mastercat[,3]*xstep*ystep}
 }
 if(peakscale){
  currentpeak = max(mastercat[, 3])
       mastercat[, 3] = mastercat[, 3] * peak/currentpeak
 }
}else{
 mastercat=cbind(x=x,y=y,mastercat=0)
 if(peakscale){
  mastercat[x+xstep/2>xcen & x-xstep/2<xcen & y+ystep/2>ycen & y-ystep/2<ycen,3]=1
 }
}
return=mastercat
}


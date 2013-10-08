inap=function(x,y,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=T){
if(majax>0){
x=x-xcen;y=y-ycen;r=sqrt(x^2+y^2)
if(deg){ang=atan2(x,y)-axang*pi/180}else{ang=atan2(x,y)-axang}
x=r*sin(ang);y=r*cos(ang)
x=x/axrat
x=x/majax;y=y/majax
return(x^2+y^2<=1)
}else{
return(x==xcen & y==ycen)
}
}

checkgrid=function(x,y,xstep,ystep,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=T){
MxMy=inap(x,y,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
LxLy=inap(x-xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
LxHy=inap(x-xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
HxLy=inap(x+xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
HxHy=inap(x+xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
full= LxLy & LxHy & HxLy & HxHy & MxMy
LxLy=inap(x-xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
LxHy=inap(x-xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
HxLy=inap(x+xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
HxHy=inap(x+xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
part= ((LxLy | LxHy | HxLy | HxHy) & full==F)
empty= (! full) & (! part)

#if(all(part==FALSE)){(part=xcen<x+xstep/2 & xcen>x-xstep/2 & ycen<y+ystep/2 & ycen>y-ystep/2)} #this is to catch situations where the aperture is smaller than the pixel grid
Npart= LxLy + LxHy + HxLy + HxHy
return=list(logic=cbind(part=part,empty=empty,full=full,mid=MxMy),Npart=Npart)
}

reexpand.grid=function(x,y,xstep,ystep,id,upres=2){
newxstep=xstep/upres
newystep=ystep/upres
tempxgrid=seq(-xstep/2+newxstep/2,xstep/2-newxstep/2,length=upres)
tempygrid=seq(-ystep/2+newystep/2,ystep/2-newystep/2,length=upres)
innergrid=expand.grid(tempxgrid,tempygrid)
loc=expand.grid(1:length(x),1:length(innergrid[,1]))
if(missing(id)){
return=list(x=x[loc[,1]]+innergrid[loc[,2],1],y=y[loc[,1]]+innergrid[loc[,2],2],xstep=xstep/upres,ystep=ystep/upres)
}else{
return=list(x=x[loc[,1]]+innergrid[loc[,2],1],y=y[loc[,1]]+innergrid[loc[,2],2],id=id[loc[,1]],xstep=xstep/upres,ystep=ystep/upres)
}
}

iterapint=function(x,y,xstep,ystep,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=T,upres=2,itersteps=9,pixscale=FALSE,peakscale=FALSE,peak=1){

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





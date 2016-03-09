reexpand.grid <-
function(x,y,xstep,ystep,id,resample.upres=2){
newxstep=xstep/resample.upres
newystep=ystep/resample.upres
tempxgrid=seq(-xstep/2+newxstep/2,xstep/2-newxstep/2,length=resample.upres)
tempygrid=seq(-ystep/2+newystep/2,ystep/2-newystep/2,length=resample.upres)
innergrid=expand.grid(tempxgrid,tempygrid)
loc=expand.grid(1:length(x),1:length(innergrid[,1]))
if(missing(id)){
return=list(x=x[loc[,1]]+innergrid[loc[,2],1],y=y[loc[,1]]+innergrid[loc[,2],2],xstep=xstep/resample.upres,ystep=ystep/resample.upres)
}else{
return=list(x=x[loc[,1]]+innergrid[loc[,2],1],y=y[loc[,1]]+innergrid[loc[,2],2],id=id[loc[,1]],xstep=xstep/resample.upres,ystep=ystep/resample.upres)
}
}

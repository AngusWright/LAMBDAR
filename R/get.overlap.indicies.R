get.overlap.indicies<-function(astr1,astr2,expand=FALSE) { 
  arr2<-(ad.to.xy(xy.to.ad(c(1.5,astr1$NAXIS[1]+0.5),c(1.5,astr1$NAXIS[2]+0.5),astr=astr1),astr=astr2))
  arr1<-(ad.to.xy(xy.to.ad(c(1.5,astr2$NAXIS[1]+0.5),c(1.5,astr2$NAXIS[2]+0.5),astr=astr2),astr=astr1))
  if (any(arr1<1)) { arr1[which(arr1<1)]<-1.5 }
  if (any(arr1[,1]>astr1$NAXIS[1])) { arr1[which(arr1[,1]>astr1$NAXIS[1]),1]<-astr1$NAXIS[1]+0.5 }
  if (any(arr1[,2]>astr1$NAXIS[2])) { arr1[which(arr1[,2]>astr1$NAXIS[2]),2]<-astr1$NAXIS[2]+0.5 }
  if (any(arr2<1)) { arr2[which(arr2<1)]<-1.5 }
  if (any(arr2[,1]>astr2$NAXIS[1])) { arr2[which(arr2[,1]>astr2$NAXIS[1]),1]<-astr2$NAXIS[1]+0.5 }
  if (any(arr2[,2]>astr2$NAXIS[2])) { arr2[which(arr2[,2]>astr2$NAXIS[2]),2]<-astr2$NAXIS[2]+0.5 }

  arr1<-floor(zapsmall(arr1))
  arr2<-floor(zapsmall(arr2))
  if (expand) { 
    XY1<-expand.grid(seq(arr1[1,1],arr1[2,1]),seq(arr1[1,2],arr1[2,2]))
    XY2<-expand.grid(seq(arr2[1,1],arr2[2,1]),seq(arr2[1,2],arr2[2,2]))
    XY1<-cbind(XY1$Var1,XY1$Var2)
    XY2<-cbind(XY2$Var1,XY2$Var2)
  } else {
    XY1<-list(X=seq(arr1[1,1],arr1[2,1]),Y=seq(arr1[1,2],arr1[2,2]))
    XY2<-list(X=seq(arr2[1,1],arr2[2,1]),Y=seq(arr2[1,2],arr2[2,2]))
  }
  return=list(arr1.XY=XY1,arr2.XY=XY2)
}


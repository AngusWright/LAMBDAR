rotate.data.2d<-function (x, y, theta,about.centre=FALSE) {
  out = make.rotation.matrix.2d(theta) %*% rbind(x, y)
  if (about.centre) { 
    out[1,]<-out[1,]+median(x)-median(out[1,])
    out[2,]<-out[2,]+median(y)-median(out[2,])
  }
  return = cbind(out[1, ], out[2, ])
}
make.rotation.matrix.2d <-function (theta) {
  theta = theta * pi/180
  sintheta = sin(theta)
  costheta = cos(theta)
  return = matrix(c(costheta, -sintheta, sintheta, costheta), ncol = 2, byrow = TRUE)
}

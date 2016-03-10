rotate.data.2d<-function (x, y, theta) {
  out = make.rotation.matrix.2d(theta) %*% rbind(x, y)
  return = cbind(out[1, ], out[2, ])
}
make.rotation.matrix.2d <-function (theta) {
  theta = theta * pi/180
  sintheta = sin(theta)
  costheta = cos(theta)
  return = matrix(c(costheta, -sintheta, sintheta, costheta), ncol = 2, byrow = TRUE)
}

interp2D <-
function(x,y,obj){
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    nx = length(xobj)
    ny = length(yobj)
    lx = approx(xobj, 1:nx, x, rule=2)$y
    ly = approx(yobj, 1:ny, y, rule=2)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    ey = ly - ly1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    return = cbind(X=x,Y=y,Z=
	zobj[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
	zobj[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
	zobj[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
	zobj[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}

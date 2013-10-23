lastnchar <-
function(x, n){
  #Return the last N characters of a string
  substr(x, nchar(x)-n+1, nchar(x))
}

 getNthVar<-
 function(varList,n,inenv,outenv,lim) {
   if (identical(inenv,outenv)) { stop("Input and Output Environments must differ") }
   for (i in varList) {
     var<-get(i, envir=inenv)
     if (length(var)>1) { var<-c(t(var)) }
     if (length(var)!=lim) { assign(i, var[1], envir=outenv) } else { assign(i, var[n], envir=outenv) }
   }
 }


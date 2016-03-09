#
# Function for returning parameters to Loop
#
getNthVar<-
 function(varList,n,inenv,outenv,lim) {
   if (identical(inenv,outenv)) { stop("Input and Output Environments must differ") }
   if (any(lsos(env=inenv)$Type=="function")) {
     warning("Function variables are not passed between environments")
   }
   #Loop through Variable list
   for (i in varList) {
     #If variable doesn't exist, set to null and loop
     if (!exists(i,envir=inenv)) {
       assign(i,NULL,envir=outenv)
       next
     }
     #If variable is a function, loop
     if (lsos(i,envir=inenv)$Type=="function") { next }
     #get variable
     var<-get(i, envir=inenv)
     if (length(var)>1) { var<-c(t(var)) }
     if (length(var)==0) {
       #Variable is NULL
       assign(i,NULL,envir=outenv)
     } else if (length(var)==lim) {
       #If variable is of full length, use the nth entry
       assign(i, var[n], envir=outenv)
     } else if (lim%%length(var)==0){
       #If variable is a multiple of the full length, use the lim%%nth entry
       ind<-n%%length(var)
       assign(i, var[ifelse(ind==0,length(var),ind)], envir=outenv)
     } else {
       #If variable isn't a multiple of full length, message a warning
       message(paste("WARNING: Variable",i,"has length that is not a factor of the longest variable (length",lim,"). Using _first_ entry."))
       assign(i, var[1], envir=outenv)
     }
   }
 }


#
#
# Function to print functions 
# that are currently flagged for debugging
# Stolen from Stackoverflow user cbeleites: 
# https://stackoverflow.com/questions/6950602/listing-functions-with-debug-flag-set-in-r
# with 1 minor addition
#
ls.deb  <- function(items = search (),simple=FALSE){
  .ls.deb <-  function (i){
    f <- ls (i)
    f <- mget (f, as.environment (i), mode = "function",

               ## return a function that is not debugged
               ifnotfound = list (function (x) function () NULL)
               )

    if (length (f) == 0)
      return (NULL)

    f <- f [sapply (f, isdebugged)]
    f <- names (f)

    ## now check whether the debugged function is masked by a not debugged one
    masked <- !sapply (f, function (f) isdebugged (get (f)))

    ## generate pretty output format:
    ## "package::function"  and "(package::function)" for masked debugged functions
    if (length (f) > 0) {
      if (!simple) { 
        if (grepl ('^package:', i)) {
          i <- gsub ('^package:', '', i)
          f <- paste (i, f, sep = "::")
        }
        f [masked] <- paste ("(", f [masked], ")", sep = "")
      }

      f
    } else {
      NULL
    }
  }
  functions <- lapply (items, .ls.deb)
  unlist (functions)
}

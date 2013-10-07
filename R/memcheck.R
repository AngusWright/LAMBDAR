object.sizes <- function(obs=ls())
{
    for (i in 1:length(obs)) {
    	paste(cat(paste(obs[i],": ")),print(object.size(obs[i]),units="auto"))
    }
}

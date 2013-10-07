make.fits.header <- function(hdr){
    
    # columns
    key = hdr[,"key"]
    value = hdr[,"value"]
    comment = hdr[,"comment"]
    
    # loop over each key
    rows = {}
    for(i in 1:length(key)){
    
        # special key?
        special = FALSE
        if(key[i]=="COMMENT" | key[i]=="HISTORY"){
            special = TRUE
        }
        
        # key
        key[i] = sprintf("%-8s", key[i])
        hier = FALSE
        if(nchar(key[i])>8){
            key[i] = paste("HIERARCH",key[i],"")
            hier = TRUE
        }
        
        # value
        if(special){
            value[i] = ""
        }else{
            # treat as a character?
            char = FALSE
            if(is.na(value[i])){
                char = TRUE
                value[i] = paste("'",sprintf("%-8s", value[i]),"'",sep="")
            }else if(is.na(suppressWarnings(as.numeric(value[i]))) & value[i]!="T" & value[i]!="F"){
                char = TRUE
                value[i] = paste("'",sprintf("%-8s", value[i]),"'",sep="")
            }
            # hierarch or normal?
            if(hier){
                pad = 33-nchar(key[i])-2-nchar(value[i])-3
                if(pad<0){pad=0}
                if(char){
                    value[i] = paste("= ", value[i],paste(rep(" ",pad),sep="",collapse=""), " / ", sep="", collapse="")
                }else{
                    value[i] = paste("= ", paste(rep(" ",pad),sep="",collapse=""),value[i], " / ", sep="", collapse="")
                }
            }else{
                if(char){
                    value[i] = paste("= ", sprintf("%-20s", value[i]), " / ", sep="")
                }else{
                    value[i] = paste("= ", sprintf("%20s", value[i]), " / ", sep="")
                }
            }
        }
        
        # update rows
        row = substr(sprintf("%-80s", paste(key[i], value[i], comment[i], sep="")),1,80)
        rows = c(rows,row)
        
    }
    
    # pad header to multiple of 2880 characters (36*80)
    hdr = paste(paste(rows,collapse="",sep=""),"END",sep="")
    pad = paste(rep(" ",(2880 - (nchar(hdr)%%2880))),collapse="",sep="")
    hdr = paste(hdr, pad, collapse="", sep="")
    
    return(hdr)
    
}


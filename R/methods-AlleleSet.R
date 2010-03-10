setMethod("bothStrands", "AlleleSet",
          function(object){
            grp1 <- c("alleleA", "alleleB")
            grp2 <- c("senseAlleleA", "senseAlleleB",
                      "antisenseAlleleA", "antisenseAlleleB")
            elem <- assayDataElementNames(object)
##            if (all(elem %in% grp1)){
	    if(all(grp1 %in% elem)){
		    return(FALSE)
            }else if (all(grp2 %in% elem)){
              return(TRUE)
            }else{
              stop("Invalid 'AlleleSet' object.")
            }
    })

setMethod("allele", "AlleleSet",
          function(object, allele, strand){
            stopifnot(!missing(allele))
            allele <- match.arg(allele, c("A", "B"))
            both <- bothStrands(object)
            if (!both){
              what <- paste("allele", allele, sep="")
            }else{
              stopifnot(!missing(strand))
              strand <- match.arg(strand, c("sense", "antisense"))
              what <- paste(strand, "Allele", allele, sep="")
            }
            assayDataElement(object, what)
          })

setMethod("getM", "AlleleSet",
          function(object){
            both <- bothStrands(object)
            if (!both){
              tmp <- allele(object, "A")-allele(object, "B")
            }else{
              tmp <- array(NA, dim=c(dim(object), 2),
                           dimnames=list(featureNames(object),
                             sampleNames(object),
                             c("antisense", "sense")))
              tmp[,,1] <- allele(object, "A", "antisense")-allele(object, "B", "antisense")
              tmp[,,2] <- allele(object, "A", "sense")-allele(object, "B", "sense")
            }
            return(tmp)
          })

setMethod("getA", "AlleleSet",
          function(object){
            both <- bothStrands(object)
            if (!both){
              tmp <- (allele(object, "A")+allele(object, "B"))/2
            }else{
              tmp <- array(NA, dim=c(dim(object), 2),
                           dimnames=list(featureNames(object),
                             sampleNames(object),
                             c("antisense", "sense")))
              tmp[,,1] <- (allele(object, "A", "antisense")+allele(object, "B", "antisense"))/2
              tmp[,,2] <- (allele(object, "A", "sense")+allele(object, "B", "sense"))/2
            }
            return(tmp)
          })

setMethod("db", "AlleleSet", function(object) db(get(annotation(object))))


setMethod("A", "AlleleSet", function(object) allele(object, "A"))
setMethod("B", "AlleleSet", function(object) allele(object, "B"))
setReplaceMethod("A", "AlleleSet", function(object, value) assayDataElementReplace(object, "A", value))
setReplaceMethod("B", "AlleleSet", function(object, value) assayDataElementReplace(object, "B", value))

setMethod("db", "AlleleSet", function(object) db(get(annotation(object))))
setMethod("A", "AlleleSet", function(object, ...) allele(object, "A", ...))
setMethod("B", "AlleleSet", function(object, ...) allele(object, "B", ...))
setReplaceMethod("A", "AlleleSet", function(object, value) {
	assayDataElementReplace(object, "alleleA", value)
})
setReplaceMethod("B", "AlleleSet", function(object, value) {
	assayDataElementReplace(object, "alleleB", value)
})

setMethod("bothStrands", "AlleleSet",
          function(object){
            grp1 <- c("alleleA", "alleleB")
            grp2 <- c("senseAlleleA", "senseAlleleB",
                      "antisenseAlleleA", "antisenseAlleleB")
            elem <- assayDataElementNames(object)
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
            ffmat <- all(unlist(eapply(assayData(object), is.ffmatrix)))
            ismat <- all(unlist(eapply(assayData(object), is.matrix)))
            stopifnot(ffmat || ismat)
            
            if (!both){
              if (ismat){
                tmp <- A(object)-B(object)
              }else{
                tmp <- ff(vmode="double", dim=dim(object))
                for (i in 1:ncol(object))
                  tmp[,i] <- A(object)[,i]-B(object)[,i]
              }
            }else{
              if (ismat){
                tmp <- array(NA, dim=c(dim(object), 2),
                             dimnames=list(featureNames(object),
                               sampleNames(object),
                               c("antisense", "sense")))
                tmp[,,1] <- A(object, "antisense")-B(object, "antisense")
                tmp[,,2] <- A(object, "sense")-B(object, "sense")
              }else{
                tmp <- ff(vmode="double", dim=c(dim(object), 2))
                for (i in 1:ncol(object)){
                  tmp[, i, 1] <- A(object, "antisense")[,i]-B(object, "antisense")[,i]
                  tmp[, i, 2] <- A(object, "sense")[,i]-B(object, "sense")[,i]
                }
              }
            }
            return(tmp)
          })

setMethod("getA", "AlleleSet",
          function(object){
            both <- bothStrands(object)
            ffmat <- all(unlist(eapply(assayData(object), is.ffmatrix)))
            ismat <- all(unlist(eapply(assayData(object), is.matrix)))
            stopifnot(ffmat || ismat)
##            rm(ffmat, ismat)
            
            if (!both){
              if (ismat){
                tmp <- (A(object)+B(object))/2
              }else{
                tmp <- ff(vmode="double", dim=dim(object))
                for (i in 1:ncol(object))
                  tmp[,i] <- (A(object)[,i]+B(object)[,i])/2
              }
            }else{
              if (ismat){
                tmp <- array(NA, dim=c(dim(object), 2),
                             dimnames=list(featureNames(object),
                               sampleNames(object),
                               c("antisense", "sense")))
                tmp[,,1] <- (A(object, "antisense")+B(object, "antisense"))/2
                tmp[,,2] <- (A(object, "sense")+B(object, "sense"))/2
              }else{
                tmp <- ff(vmode="double", dim=c(dim(object), 2))
                for (i in 1:ncol(object)){
                  tmp[, i, 1] <- (A(object, "antisense")[,i]+B(object, "antisense")[,i])/2
                  tmp[, i, 2] <- (A(object, "sense")[,i]+B(object, "sense")[,i])/2
                }
              }
            }
            return(tmp)
          })





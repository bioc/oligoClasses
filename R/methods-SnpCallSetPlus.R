##setMethod("getA", "SnpCallSetPlus",
##          function(obj){
##            tmp <- array(NA, dim=c(nrow(antisenseThetaA(obj)),
##                               ncol(antisenseThetaA(obj)), 2),
##                         dimnames=list(rownames(antisenseThetaA(obj)),
##                           colnames(antisenseThetaA(obj)),
##                           c("antisense", "sense")))
##            tmp[,,1] <- .5*(antisenseThetaA(obj)+antisenseThetaB(obj))
##            tmp[,,2] <- .5*(senseThetaA(obj)+senseThetaB(obj))
##            return(tmp)
##          })


setMethod("db", "SnpCallSetPlus",
          function(object) {
		  require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })




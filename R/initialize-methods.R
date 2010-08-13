setMethod("initialize", "CopyNumberSet",
          function(.Object,
                   assayData = assayDataNew(copyNumber = copyNumber,
                                            cnConfidence = cnConfidence, ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   protocolData = phenoData[,integer(0)],
                   copyNumber = new("matrix"),
                   cnConfidence = matrix(numeric(),
                                            nrow=nrow(copyNumber), ncol=ncol(copyNumber),
                                            dimnames=dimnames(copyNumber)),
                   ...) {
		  .Object <- callNextMethod(.Object,
					    assayData = assayData,
					    phenoData = phenoData,
					    featureData = featureData,
					    experimentData = experimentData,
					    annotation = annotation,
					    protocolData = protocolData)
		  if(checkAnnotation(annotation))
			  .Object <- annotate(.Object)
		  return(.Object)
          })

setMethod("initialize", "oligoSnpSet",
	  function(.Object,
		   call=new("matrix"),
		   callProbability=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   copyNumber=matrix(numeric(), nrow=nrow(call), ncol=ncol(call),  dimnames=dimnames(call)),
		   cnConfidence=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   annotation=character(), ...){
		  .Object <- callNextMethod(.Object,
					    call=call,
					    callProbability=callProbability,
					    copyNumber=copyNumber,
					    cnConfidence=cnConfidence,
					    annotation=annotation, ... )
		  if(checkAnnotation(annotation))
			  .Object <- annotate(.Object)  		  
		  return(.Object)
	  })

setValidity("oligoSnpSet", function(object) {
	assayDataValidMembers(assayData(object), c("call", "callProbability", "copyNumber", "cnConfidence"))
})

## RS: ask BC about this... initialization method for CNSet does not work when this is uncommented
setValidity("AlleleSet",
            function(object){
              grp1 <- c("alleleA", "alleleB")
              grp2 <- c("senseAlleleA", "senseAlleleB",
                        "antisenseAlleleA", "antisenseAlleleB")
              elem <- assayDataElementNames(object)
              ok <- all(grp1 %in% elem) || all(grp2 %in% elem)
              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
              if (!ok){
                paste("Elements of 'AlleleSet' must be:",
                      f(grp1), "OR", f(grp2))
              }else{
                TRUE
              }
            })

##setMethod("initialize", "AlleleSet",
##	  function(.Object, ...){
##		  args <- names(list(...))
##		  grp1 <- c("alleleA", "alleleB")
##		  grp2 <- c("senseAlleleA", "senseAlleleB",
##			    "antisenseAlleleA", "antisenseAlleleB")
##		  f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
##		  if(all(grp1 %in% args)){
##			  browser()
##		  }
##		  if(all(grp2 %in% args)){
##			  
##		  }
##	  })
##		  browser()
##		  if(any(args%in%grp1)){
##			  browser()
##			  .Object <- callNextMethod(.Object, ...)
##		  } else{
##			  stop(paste("Elements of 'AlleleSet' must be:",
##				f(grp1), "OR", f(grp2)))
##		  }
##		  return(.Object)
##	  })

## BC: AlleleSet must allow for sense/antisense chips too
## BC: I'm commenting out this initialization method
## BC: to allow the standard eSet initialization
## BC: therefore, all new instances must be named:
## BC: eg: new("AlleleSet", alleleA=<value>, alleleB=<value>)
## setMethod("initialize", "AlleleSet",
##           function(.Object,
## 		   alleleA=new("matrix"),
## 		   alleleB=new("matrix"), ...){
## 		  callNextMethod(.Object, alleleA=alleleA, alleleB=alleleB, ...)
## 	  })

setMethod("initialize", "SnpSuperSet",
          function(.Object, call=new("matrix"), callProbability=new("matrix"), ...){
		  callNextMethod(.Object, call=call, callProbability=callProbability, ...)
	  })

## b/c I removed alleleA/alleleB from AlleleSet initialization, I'll add it here
## as this seems to be the place where such arguments are expected.
setMethod("initialize", "CNSet",
          function(.Object,
		   ##CA=new("matrix"),
		   ##CB=new("matrix"),
		   alleleA=new("matrix"),
		   alleleB=new("matrix"),
##		   batch=as.factor(rep(1, ncol(CA))),
		   ...
		   ){
		  .Object <- callNextMethod(.Object,
					    ##CA=CA,
					    ##CB=CB,
					    alleleA=alleleA,
					    alleleB=alleleB,
					    ...)
##		  if(length(batch) == ncol(.Object)){
##			  .Object@batch <- batch
##		  } else stop("batch variable must be same length as number of samples")
          })

setValidity("CNSet", function(object) {
	assayDataValidMembers(assayData(object), c("call", "callProbability", "alleleA", "alleleB"))
})

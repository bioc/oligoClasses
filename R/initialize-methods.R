setMethod("initialize", "oligoSnpSet",
	  function(.Object, assayData,
		   call=new("matrix"),
		   callProbability=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   copyNumber=matrix(numeric(), nrow=nrow(call), ncol=ncol(call),  dimnames=dimnames(call)),
		   cnConfidence=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)), ... ){
		  if(missing(annotation)){
			  stop("must specify annotation")
		  } else{
			  stopifnot(isValidAnnotation(annotation))
			  .Object@annotation <- annotation
		  }
		  if(missing(assayData)){
			  .Object <- callNextMethod(.Object,
						    call=call,
						    callProbability=callProbability,
						    copyNumber=copyNumber,
						    cnConfidence=cnConfidence, ...)
		  } else{
			  .Object <- callNextMethod(.Object,
						    assayData=assayData, ...)
		  }
		  ## Do after annotation has been assigned
		  if(!(all(c("chromosome", "position", "isSnp") %in% fvarLabels(.Object)))){
			  .Object@featureData <- addFeatureAnnotation.eSet(.Object)
		  }
	  })

isValidAnnotation <- function(x){
	validAnn <- validAnnotation()
	x %in% validAnn
}

validAnnotation <- function(){
	c("pd.mapping50k.hind240", "pd.mapping50k.xba240",
	  "pd.mapping50k.hind240,pd.mapping50k.xba240",	  
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp,pd.mapping250k.sty",
	  "pd.genomewidesnp.5",
	  "pd.genomewidesnp.6")
}


setValidity("oligoSnpSet", function(object) {
	assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "copyNumber", "cnConfidence"))
})

setValidity("AlleleSet",
            function(object){
              grp1 <- c("alleleA", "alleleB")
              grp2 <- c("senseAlleleA", "senseAlleleB",
                        "antisenseAlleleA", "antisenseAlleleB")
              elem <- assayDataElementNames(object)
              ok <- all(grp1 %in% elem) || all(grp2 %in% elem)
              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
              if (!ok){
                paste("Elements of 'AlleleSummarySet' must be:",
                      f(grp1), "OR", f(grp2))
              }else{
                TRUE
              }
            })

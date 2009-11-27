setMethod("initialize", "oligoSnpSet",
	  function(.Object, assayData, calls=new("matrix"),
		   callsConfidence=matrix(numeric(), nrow=nrow(calls), ncol=ncol(calls), dimnames=dimnames(calls)),
		   copyNumber=matrix(numeric(), nrow=nrow(calls), ncol=ncol(calls),  dimnames=dimnames(calls)),
		   cnConfidence=matrix(numeric(), nrow=nrow(calls), ncol=ncol(calls), dimnames=dimnames(calls)), ... ){
	if(missing(assayData)){
		callNextMethod(.Object,
			       calls=calls,
			       callsConfidence=callsConfidence,
			       copyNumber=copyNumber,
			       cnConfidence=cnConfidence, ...)
	} else{
		callNextMethod(.Object,
			       assayData=assayData, ...)
	}
})

setValidity("oligoSnpSet", function(object) {
	assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "copyNumber", "cnConfidence"))
})

setValidity("AlleleSet",
            function(object){
              grp1 <- c("alleleA", "alleleB")
              grp2 <- c("senseAlleleA", "senseAlleleB",
                        "antisenseAlleleA", "antisenseAlleleB")
              elem <- assayDataElementNames(object)
              ok <- all(elem %in% grp1) || all(elem %in% grp2)
              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
              if (!ok){
                paste("Elements of 'AlleleSummarySet' must be:",
                      f(grp1), "OR", f(grp2))
              }else{
                TRUE
              }
            })

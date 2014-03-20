setMethod("batchNames", "AssayData",
	  function(object){
		  ##should call method for AssayData
		  sampleNames(object)
	  })

setReplaceMethod("batchNames", "AssayData",
		 function(object, value){
			 sampleNames(object) <- value
			 return(object)
		 })

setMethod("nu", c("AssayData", "character"),
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 A="nuA",
				 B="nuB",
				 stop("allele must be 'A' or 'B'"))
		  }
		  val <- getValue(allele)
		  assayDataElement(object, val)
})



setMethod("phi", c("AssayData", "character"),
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 A="phiA",
				 B="phiB",
				 stop("allele must be 'A' or 'B'"))
		  }
		  val <- getValue(allele)
		  assayDataElement(object, val)
	  })

##setMethod("sigma2", c("AssayData", "character"),
##	  function(object, allele){
##		  getValue <- function(allele){
##			  switch(allele,
##				 A="sig2A",
##				 B="sig2B",
##				 stop("allele must be 'A' or 'B'"))
##		  }
##		  val <- getValue(allele)
##		  assayDataElement(object, val)
##	  })

##setMethod("tau2", c("AssayData", "character"),
##	  function(object, allele){
##		  getValue <- function(allele){
##			  switch(allele,
##				 A="tau2A",
##				 B="tau2B",
##				 stop("allele must be 'A' or 'B'"))
##		  }
##		  val <- getValue(allele)
##		  assayDataElement(object, val)
##	  })

##setMethod("corr", c("AssayData", "character"),
##	  function(object, allele){
##		  getValue <- function(allele){
##			  switch(allele,
##				 AA="corrAA",
##				 AB="corrAB",
##				 BB="corrBB",
##				 stop("allele must be 'AA', 'AB', or 'BB'"))
##		  }
##		  val <- getValue(allele)
##		  assayDataElement(object, val)
##	  })

setMethod("flags", signature(object="AssayData"), function(object) assayDataElement(object, "flags"))

AssayDataList <- function(storage.mode = c("lockedEnvironment", "environment", "list"), ...) {
  storage.mode <- match.arg(storage.mode) ## defaults to "lockedEnvironment"
  assayData <- switch(storage.mode,
                      lockedEnvironment =,
                      environment = new.env(parent=emptyenv()),
                      list = list())
  arglist <- list(...)
  for (nm in names(arglist)) assayData[[nm]] <- arglist[[nm]]
  ## FIX:  ::: not safe
  if (storage.mode == "lockedEnvironment") Biobase:::assayDataEnvLock(assayData)
  assayData
}

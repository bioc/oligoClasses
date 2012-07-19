## utilities for parallel computing (now (20Mar2012) via foreach)
##
## NOTE: Using this framework, it is possible to allow for parallel
## computing WITHOUT the ff package. The problem with this is to
## overload the system as operations will pottentially copy the same
## data over and over and over. Therefore, I chose to allow parallel
## computing only when ff is loaded.
##
## Summary (useful when coding):
##   - parStatus: TRUE if requirements for parallel are met
##   - ocProbesets: number of probesets to process at a time (batch)
##   - ocSamples: number of samples to process at a time (batch)
##   - ocPath: path where ff objects are to be saved

parStatus <- function()
    getDoParRegistered() & isPackageLoaded('ff')

ocParallelStatus <- function(verbose=TRUE){
  cl <- parStatus()
  if (verbose){
      message("Parallel computing support for 'oligo/crlmm': ", appendLF=FALSE)
      if (!cl){
          message("Disabled")
          message("     - Load 'ff'")
          message("     - Load and register a 'foreach' adaptor")
          message("        Example - Using 'multicore' for 2 cores:")
          message("             library(doMC)")
          message("             registerDoMC(2)")
      }else{
          message("Enabled")
          ocProbesets(getOption('ocProbesets'))
          ocSamples(getOption('ocSamples'))
      }
      message(getBar())
  }
  return(cl)
}

ocProbesets <- function(n){
  if (missing(n)){
    return(getOption("ocProbesets"))
  }else{
    options(ocProbesets=n)
    invisible(TRUE)
  }
}

ocSamples <- function(n){
  if (missing(n)){
    return(getOption("ocSamples"))
  }else{
    options(ocSamples=n)
    invisible(TRUE)
  }
}

ocLapply <- function(X, FUN, ..., neededPkgs){
    if(missing(neededPkgs)) neededPkgs <- 'ff'
    else neededPkgs <- unique(c('ff', neededPkgs))
    x <- NULL ## to make NOTE go away in R's package checker
    if (parStatus()){
       res <- foreach(x=X, .packages=neededPkgs) %dopar% FUN(x, ...)
    }else{
       res <- lapply(X, FUN, ...)
    }
  return(res)
}

splitIndicesByLength <- function(x, lg, balance=FALSE){
	lx <- length(x)
	split(x, rep(seq(1,lx), each=lg, length.out=lx))
}

splitIndicesByNode <- function(x){
    split(x, sort(rep(1:getDoParWorkers(), length.out=length(x))))
}

## deprecated
setCluster <- function(...){
    msg <- paste('To set cluster for oligo/crlmm/friends,',
                 'use the "foreach" package and one adaptor like "doMC" or "doMPI."',
                 'Then, register the adaptor via registerDo*().')
    .Deprecated(msg=msg)
}

delCluster <- function(){
    msg <- paste('Use the recommendations used by the "foreach" package.')
    .Deprecated(msg=msg)
}


getCluster <- function(){
    .Deprecated('getDoParWorkers')
}

requireClusterPkgSet <- function(packages){
    .Deprecated(msg='Function no longer needed. Replaced by proper argument of foreach')
}

requireClusterPkg <- function(pkg, character.only=TRUE){
    .Deprecated(msg='Function no longer needed. Replaced by proper argument of foreach')
}

## utilities for parallel computing (via snow)
##
## Summary (useful when coding):
##   - parStatus: TRUE if requirements for parallel are met
##   - ocProbesets: number of probesets to process at a time (batch)
##   - ocSamples: number of samples to process at a time (batch)
##   - ocPath: path where ff objects are to be saved

parStatus <- function()
    is(getOption("cluster"), "cluster") && isPackageLoaded("snow")

ocParallelStatus <- function(verbose=TRUE){
  sn <- isPackageLoaded("snow")
  cl <- parStatus()
  ld <- isPackageLoaded("ff")
  if (verbose){
    message("Parallel computing support for 'oligo/crlmm': ", appendLF=FALSE)
    if (!ld){  
      message("Disabled")
      message("     - Load 'ff'")
      if (!sn){
        message("     - Load 'snow'")
        message("     - Use options(cluster=makeCluster(...))")
      } else {
        if (!cl)
          message("     - Use options(cluster=makeCluster(...))")
      }
    }else{
      if (sn){
        if (cl){
          message("Enabled")
        }else{
          message("Disabled")
          message("     - Use options(cluster=makeCluster(...))")
        }
      }else{
        message("Disabled")
        message("     - Load 'snow'")
        message("     - Use options(cluster=makeCluster(...))")
      }
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

setCluster <- function(...){
  pkg <- "snow"
  require(pkg, character.only=TRUE)
  options(cluster=makeCluster(...))
}

delCluster <- function(){
  stopCluster(getOption("cluster"))
  options(cluster=NULL)
}

getCluster <- function()
  getOption("cluster")

requireClusterPkgSet <- function(packages){
  if (!parStatus())
    stop("cluster is not ready. Use 'setCluster'.")
  for (pkg in packages){
    pkgOnCluster <- requireClusterPkg(pkg, character.only=TRUE)
    if (!pkgOnCluster){
      msg <- paste("Package '", pkg, "' not found on the cluster. ",
                   "Install it or load it manually using ",
                   "'clusterEvalQ(getCluster(), library(", pkg,
                   ", lib.loc=<APPROPRIATE PATH>))'", sep="")
      stop(msg)
    }
  }
  TRUE
}

requireClusterPkg <- function(pkg, character.only=TRUE)
  all(unlist(clusterCall(getCluster(), require, pkg, character.only=character.only)))

ocLapply <- function(X, FUN, ..., neededPkgs){
  if (parStatus()){
    if (missing(neededPkgs)){
      neededPkgs <- "ff"
    }else{
      neededPkgs <- unique(append(neededPkgs, "ff"))
    }
    ok <- requireClusterPkgSet(neededPkgs)
    res <- parLapply(getCluster(), X, FUN, ...)
  }else{
    res <- lapply(X, FUN, ...)
  }
  return(res)
}

splitIndicesByLength <- function(x, lg){
  lx <- length(x)
  split(x, rep(1:lx, each=lg, length.out=lx))
}

splitIndicesByNode <- function(x){
  if (parStatus()){
    clusterSplit(getCluster(), x)
  }else{
    list(x)
  }
}


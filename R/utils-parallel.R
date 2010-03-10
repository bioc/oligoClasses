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
        message("     - Use options(cluster=makeCluster(...)")
      } else {
        if (!cl)
          message("     - Use options(cluster=makeCluster(...)")
      }
    }else{
      if (sn){
        if (cl){
          message("Enabled")
        }else{
          message("Disabled")
          message("     - Use options(cluster=makeCluster(...)")
        }
      }else{
        message("Disabled")
        message("     - Load 'snow'")
        message("     - Use options(cluster=makeCluster(...)")
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

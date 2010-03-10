## utilities for Large Dataset Support
##
## Summary:
##   - is.ffmatrix: test specifically for ff_matrix
##   - ldStatus: TRUE if Large Dataset Support is available
##   - ldPath: where to save ff files
##   - createFF: creates an ff object setting path appropriately
##               (leaving out of fftempdir b/c parallel processes
##               access the object very easily)

createFF <- function(name, dim, vmode="double")
  ff(vmode=vmode, dim=dim, pattern=file.path(ldPath(), basename(name)))

setMethod("annotatedDataFrameFrom", "ff_matrix",
          Biobase:::annotatedDataFrameFromMatrix)

is.ffmatrix <- function(object)
  is(object, "ff_matrix")

ldPath <- function(path){
  if (missing(path)){
    return(getOption("ldPath"))
  }else{
    stopifnot(is.character(path))
    options(ldPath=path)
  }
}

ldSetOptions <- function(nsamples=100, nprobesets=1000,
                         path=getwd(), verbose=FALSE){
  ocProbesets(nprobesets)
  ocSamples(nsamples)
  ldPath(path)
  ldStatus(verbose)
  TRUE
}

ldStatus <- function(verbose=FALSE){
  ld <- isPackageLoaded("ff")
  if (verbose){
    message(getBar())
    message("Large dataset support for 'oligo/crlmm': ", appendLF=FALSE)
    if (ld){
      message("Enabled")
      ns <- prettyNum(c(ocProbesets(), ocSamples()), big.mark=",")
      message("    - Probesets: ", ns[1])
      message("    - Samples..: ", ns[2])
      message("    - Path.....: ", ldPath())
    }else{
      message("Disabled")
      message("     - Load 'ff'")
    }      
    message(getBar())
  }
  return(ld)
}


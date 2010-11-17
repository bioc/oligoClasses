THISPKG <- "oligoClasses"
.oligoClassesPkgEnv <- new.env(parent=emptyenv())


.onAttach <- function(libname, pkgname) {
  version <- packageDescription("oligoClasses", field="Version")
  message(getBar())
  message("Welcome to oligoClasses version ", version)
  ldSetOptions()
  bm <- ldStatus(TRUE)
  snow <- ocParallelStatus()

  setHook(packageEvent("ff", "attach"),
          function(...){
            ldSetOptions(verbose=FALSE)
            ldStatus(TRUE)
          })

  setHook(packageEvent("snow", "attach"),
          function(...){
            ocParallelStatus()
          })

}

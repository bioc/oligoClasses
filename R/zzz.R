# Loading required libraries
THISPKG <- "oligoClasses"

.onLoad <- function(libname, pkgname) {
  require("methods")
}

.onAttach <- function(libname, pkgname) {
	message("Welcome to oligoClasses version ", packageDescription(THISPKG, field="Version"))
}

.onUnload <- function( libpath ){
	library.dynam.unload(THISPKG, libpath)
}

.oligoClassesPkgEnv <- new.env(parent=emptyenv())

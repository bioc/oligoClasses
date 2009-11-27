setMethod("show", "DBPDInfo", function(object) {
    cat("Class........:", class(object), "\n")
    cat("Manufacturer.:", manufacturer(object), "\n")
    cat("Genome Build.:", genomeBuild(object), "\n")
    cat("Chip Geometry:", geometry(object)[1], "rows x ", geometry(object)[2], "columns\n")
})

setMethod("show", "FeatureSet", function(object){
  callNextMethod()
  requireAnnotation(annotation(object))
})

setMethod("initialize", "SnpCallSet",
          function(.Object,
                   calls=new("matrix"),
                   callsConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)), ... ){
            callNextMethod(.Object,
                           calls=calls,
                           callsConfidence=callsConfidence, ...)
          })

setMethod("initialize", "SnpCopyNumberSet",
          function(.Object,
                   copyNumber=new("matrix"),
                   cnConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)), ... ){
            callNextMethod(.Object,
                           copyNumber=copyNumber,
                           cnConfidence=cnConfidence, ...)
          })

setMethod("initialize", "oligoSnpSet",
          function(.Object,
                   calls=new("matrix"),
                   callsConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   copyNumber=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)),
                   cnConfidence=matrix(numeric(), nrow=nrow(calls),
                     ncol=ncol(calls),
                     dimnames=dimnames(calls)), ... ){
            callNextMethod(.Object,
                           calls=calls,
                           callsConfidence=callsConfidence,
                           copyNumber=copyNumber,
                           cnConfidence=cnConfidence, ...)
          })

setValidity("SnpCallSet", function(object) {
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence"))
})

setValidity("SnpCopyNumberSet", function(object) {
  assayDataValidMembers(assayData(object), c("copyNumber", "cnConfidence"))
})

setValidity("oligoSnpSet", function(object) {
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "copyNumber", "cnConfidence"))
})

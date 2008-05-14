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

setMethod("initialize", "SnpCallSetPlus",
          function(.Object,
                   phenoData, featureData,
                   calls=new("matrix"),
                   callsConfidence=new("matrix"),
                   antisenseThetaA=new("matrix"),
                   antisenseThetaB=new("matrix"),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"), ... ){
            ad <- assayDataNew("lockedEnvironment",
                               calls=calls,
                               callsConfidence=callsConfidence,
                               antisenseThetaA=antisenseThetaA,
                               antisenseThetaB=antisenseThetaB,
                               senseThetaA=senseThetaA,
                               senseThetaB=senseThetaB)
            
            assayData(.Object) <- ad
            if (missing(phenoData))
              phenoData(.Object) <- annotatedDataFrameFrom(calls, byrow=FALSE)
            if (missing(featureData))
              featureData(.Object) <- annotatedDataFrameFrom(calls, byrow=TRUE)
            .Object
          })

setMethod("initialize", "SnpCnvCallSetPlus",
          function(.Object,
                   phenoData, featureData,
                   calls=new("matrix"),
                   callsConfidence=new("matrix"),
                   antisenseThetaA,
                   antisenseThetaB,
                   senseThetaA,
                   senseThetaB,
                   thetaA,
                   thetaB, ... ){
            ad <- assayDataNew("lockedEnvironment",
                               calls=calls,
                               callsConfidence=callsConfidence,
                               thetaA=thetaA,
                               thetaB=thetaB)
            assayData(.Object) <- ad
            if (missing(phenoData))
              phenoData(.Object) <- annotatedDataFrameFrom(calls, byrow=FALSE)
            if (missing(featureData))
              featureData(.Object) <- annotatedDataFrameFrom(calls, byrow=TRUE)
            .Object
          })



setMethod("initialize", "SnpCopyNumberSet",
          function(.Object,
                   copyNumber=new("matrix"),
                   cnConfidence=matrix(numeric(), nrow=nrow(copyNumber),
                     ncol=ncol(copyNumber),
                     dimnames=dimnames(copyNumber)), ... ){
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

## SnpQSet ## From oligo

setMethod("initialize", "SnpQSet",
          function(.Object,
                   assayData = assayDataNew(senseThetaA=senseThetaA,
                     senseThetaB=senseThetaB,
                     antisenseThetaA=antisenseThetaA,
                     antisenseThetaB=antisenseThetaB),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"),
                   antisenseThetaA=new("matrix"),
                   antisenseThetaB=new("matrix"),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=new("character")){
            .Object <- callNextMethod(.Object,
                                  assayData = assayDataNew(
                                    senseThetaA=senseThetaA,
                                    senseThetaB=senseThetaB,
                                    antisenseThetaA=antisenseThetaA,
                                    antisenseThetaB=antisenseThetaB),
                                  phenoData=phenoData,
                                  experimentData=experimentData,
                                  annotation=annotation)
            .Object
          })

setValidity("SnpQSet",
            function(object)
            assayDataValidMembers(assayData(object),
                                  c("senseThetaA",
                                    "senseThetaB",
                                    "antisenseThetaA",
                                    "antisenseThetaB"))
            )


## TilingQSet

setMethod("initialize", "TilingQSet",
          function(.Object,
                   assayData = assayDataNew(M=M),
                   M=new("matrix"),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentalData=new("MIAME"),
                   annotation=new("character")){
            .Object <- callNextMethod(.Object,
                                      assayData = assayDataNew(M=M),
                                      phenoData=phenoData,
                                      featureData=featureData,
                                      experimentalData=experimentalData,
                                      annotation=annotation)
            .Object
          })

setValidity("TilingQSet", function(object) assayDataValidMembers(assayData(object), "M"))


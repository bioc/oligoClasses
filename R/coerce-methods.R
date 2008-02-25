setAs("SnpCallSetPlus", "oligoSnpSet",
      function(from, to){
        cn <- calculateCopyNumber(from)
        cnConf <- matrix(NA, nrow=nrow(from), ncol=ncol(from))
        rownames(cnConf) <- featureNames(from)
        colnames(cnConf) <- sampleNames(from)
        new("oligoSnpSet",
            copyNumber=cn,
            cnConfidence=cnConf,
            calls=calls(from),
            callsConfidence=callsConfidence(from),
            experimentData=experimentData(from),
            featureData=featureData(from),
            phenoData=phenoData(from),
            annotation=annotation(from))
##            notes="copy number defined as the average abundance of the A and B alleles")
      })
  

setAs("SnpCallSetPlus", "oligoSnpSet",
      function(from, to){

	      cn <- calculateCopyNumber(from)
	      ##use inverse of across sample standard deviation as a confidence score
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
      })
  
